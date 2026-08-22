#include "collision.hpp"

#include <algorithm>
#include <cmath>

namespace rbe {

AABB computeAABB(const Shape& s, const Transform& t) {
    if (s.type == ShapeType::Sphere) {
        Vec3 r{s.radius, s.radius, s.radius};
        return {t.p - r, t.p + r};
    }
    // World extent along each world axis is |R| * halfExtents.
    Vec3 e{dot(vabs(t.R.r[0]), s.halfExtents), dot(vabs(t.R.r[1]), s.halfExtents),
           dot(vabs(t.R.r[2]), s.halfExtents)};
    return {t.p - e, t.p + e};
}

namespace {

bool sphereSphere(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb, float margin,
                  ContactManifold& m) {
    Vec3 d = tb.p - ta.p;
    float dist = length(d);
    float sep = dist - a.radius - b.radius;
    if (sep > margin) return false;
    m.normal = dist > 1e-6f ? d / dist : Vec3{0, 1, 0};
    Vec3 onA = ta.p + m.normal * a.radius;
    Vec3 onB = tb.p - m.normal * b.radius;
    m.count = 1;
    m.points[0] = {(onA + onB) * 0.5f, -sep, 0};
    return true;
}

// Box A against sphere B; normal points from the box to the sphere.
bool boxSphere(const Shape& box, const Transform& tb, const Shape& sph, const Transform& ts, float margin,
               ContactManifold& m) {
    const Vec3 e = box.halfExtents;
    Vec3 local = transpose(tb.R) * (ts.p - tb.p);
    Vec3 clamped{std::clamp(local.x, -e.x, e.x), std::clamp(local.y, -e.y, e.y), std::clamp(local.z, -e.z, e.z)};
    Vec3 delta = local - clamped;
    float distSq = lengthSq(delta);

    Vec3 nLocal, surface;
    float depth;
    if (distSq > 1e-12f) {
        float dist = std::sqrt(distSq);
        if (dist - sph.radius > margin) return false;
        nLocal = delta / dist;
        depth = sph.radius - dist;
        surface = clamped;
    } else {
        // Centre is inside the box: push out through the nearest face.
        int axis = 0;
        float best = e.x - std::fabs(local.x);
        for (int k = 1; k < 3; ++k) {
            float d = e[k] - std::fabs(local[k]);
            if (d < best) {
                best = d;
                axis = k;
            }
        }
        nLocal = {0, 0, 0};
        nLocal[axis] = local[axis] >= 0 ? 1.0f : -1.0f;
        depth = sph.radius + best;
        surface = local;
        surface[axis] = nLocal[axis] * e[axis];
    }
    m.normal = tb.R * nLocal;
    Vec3 onBox = tb.p + tb.R * surface;
    Vec3 onSphere = ts.p - m.normal * sph.radius;
    m.count = 1;
    m.points[0] = {(onBox + onSphere) * 0.5f, depth, 0};
    return true;
}

struct ClipVertex {
    Vec3 p;
    uint32_t id;
};

// Keeps the part of the polygon where dot(n, p) <= offset. Returns the new vertex count.
int clipPolygon(const ClipVertex* in, int n, ClipVertex* out, const Vec3& normal, float offset, uint32_t planeId) {
    int count = 0;
    for (int i = 0; i < n; ++i) {
        const ClipVertex& a = in[i];
        const ClipVertex& b = in[(i + 1) % n];
        float da = dot(normal, a.p) - offset;
        float db = dot(normal, b.p) - offset;
        if (da <= 0) out[count++] = a;
        if ((da < 0) != (db < 0) && da != db) {
            float t = da / (da - db);
            // New vertex: id built from the edge it sits on and the clipping plane.
            out[count++] = {a.p + (b.p - a.p) * t, 0x80u | (planeId << 4) | (a.id & 0xFu)};
        }
    }
    return count;
}

// Picks up to 4 points that keep the manifold's area: deepest, farthest from it, then the two that
// enlarge the area the most.
int reduceToFour(ManifoldPoint* pts, int n, const Vec3& normal) {
    if (n <= 4) return n;
    ManifoldPoint keep[4];
    int i0 = 0;
    for (int i = 1; i < n; ++i)
        if (pts[i].depth > pts[i0].depth) i0 = i;
    keep[0] = pts[i0];

    int i1 = -1;
    float best = -1;
    for (int i = 0; i < n; ++i) {
        float d = lengthSq(pts[i].position - keep[0].position);
        if (d > best) {
            best = d;
            i1 = i;
        }
    }
    keep[1] = pts[i1];

    int i2 = -1;
    best = -1;
    for (int i = 0; i < n; ++i) {
        float area = std::fabs(dot(cross(keep[0].position - pts[i].position, keep[1].position - pts[i].position), normal));
        if (area > best) {
            best = area;
            i2 = i;
        }
    }
    keep[2] = pts[i2];

    // Fourth point: the one on the opposite side of the triangle's edges with the largest area.
    int i3 = -1;
    best = -1;
    for (int i = 0; i < n; ++i) {
        if (i == i0 || i == i1 || i == i2) continue;
        const Vec3& p = pts[i].position;
        float a01 = dot(cross(keep[0].position - p, keep[1].position - p), normal);
        float a12 = dot(cross(keep[1].position - p, keep[2].position - p), normal);
        float a20 = dot(cross(keep[2].position - p, keep[0].position - p), normal);
        float area = std::max({-a01, -a12, -a20, 0.0f});
        if (area > best) {
            best = area;
            i3 = i;
        }
    }
    int count = 3;
    if (i3 >= 0) keep[count++] = pts[i3];
    for (int i = 0; i < count; ++i) pts[i] = keep[i];
    return count;
}

struct BoxFrame {
    Vec3 c;
    Vec3 axis[3];
    Vec3 e;
};

BoxFrame frameOf(const Shape& s, const Transform& t) {
    return {t.p, {t.R.col(0), t.R.col(1), t.R.col(2)}, s.halfExtents};
}

float projectRadius(const BoxFrame& b, const Vec3& L) {
    return b.e.x * std::fabs(dot(b.axis[0], L)) + b.e.y * std::fabs(dot(b.axis[1], L)) +
           b.e.z * std::fabs(dot(b.axis[2], L));
}

// Face contact with `ref` as the reference box. nRef points from ref toward inc.
void faceContact(const BoxFrame& ref, int refAxis, const Vec3& nRef, const BoxFrame& inc, float margin,
                 uint32_t idBase, ContactManifold& m, bool flip) {
    // Incident face: the face of `inc` most anti-parallel to nRef.
    int incAxis = 0;
    float best = -1;
    for (int k = 0; k < 3; ++k) {
        float d = std::fabs(dot(inc.axis[k], nRef));
        if (d > best) {
            best = d;
            incAxis = k;
        }
    }
    Vec3 incNormal = dot(inc.axis[incAxis], nRef) > 0 ? -inc.axis[incAxis] : inc.axis[incAxis];
    Vec3 incCenter = inc.c + incNormal * inc.e[incAxis];
    int k1 = (incAxis + 1) % 3, k2 = (incAxis + 2) % 3;
    Vec3 u = inc.axis[k1] * inc.e[k1], v = inc.axis[k2] * inc.e[k2];

    ClipVertex poly[16], tmp[16];
    poly[0] = {incCenter + u + v, 0};
    poly[1] = {incCenter - u + v, 1};
    poly[2] = {incCenter - u - v, 2};
    poly[3] = {incCenter + u - v, 3};
    int n = 4;

    // Clip against the 4 side planes of the reference face.
    int r1 = (refAxis + 1) % 3, r2 = (refAxis + 2) % 3;
    const Vec3 sideAxes[4] = {ref.axis[r1], -ref.axis[r1], ref.axis[r2], -ref.axis[r2]};
    const float sideExt[4] = {ref.e[r1], ref.e[r1], ref.e[r2], ref.e[r2]};
    for (int s = 0; s < 4 && n > 0; ++s) {
        float offset = dot(sideAxes[s], ref.c) + sideExt[s];
        n = clipPolygon(poly, n, tmp, sideAxes[s], offset, static_cast<uint32_t>(s));
        std::copy(tmp, tmp + n, poly);
    }

    Vec3 refFaceCenter = ref.c + nRef * ref.e[refAxis];
    uint32_t refFace = static_cast<uint32_t>(refAxis * 2 + (dot(ref.axis[refAxis], nRef) > 0 ? 0 : 1));
    uint32_t incFace = static_cast<uint32_t>(incAxis * 2 + (dot(inc.axis[incAxis], incNormal) > 0 ? 0 : 1));

    ManifoldPoint pts[16];
    int count = 0;
    for (int i = 0; i < n; ++i) {
        float depth = -dot(nRef, poly[i].p - refFaceCenter);
        if (depth < -margin) continue;
        pts[count++] = {poly[i].p + nRef * (0.5f * depth), depth,
                        idBase | (refFace << 8) | (incFace << 12) | (flip ? 0x10000u : 0u) | poly[i].id};
    }
    m.normal = flip ? -nRef : nRef;
    if (count > 8) count = reduceToFour(pts, count, m.normal);
    m.count = count;
    for (int i = 0; i < count; ++i) m.points[i] = pts[i];
}

void closestPointsSegments(const Vec3& p1, const Vec3& d1, float h1, const Vec3& p2, const Vec3& d2, float h2,
                           Vec3& c1, Vec3& c2) {
    // Lines p1 + s*d1 and p2 + t*d2 with unit directions, s in [-h1,h1], t in [-h2,h2].
    Vec3 r = p1 - p2;
    float b = dot(d1, d2), c = dot(d1, r), f = dot(d2, r);
    float denom = 1.0f - b * b;
    float s = denom > 1e-6f ? std::clamp((b * f - c) / denom, -h1, h1) : 0.0f;
    float t = std::clamp(b * s + f, -h2, h2);
    s = std::clamp(b * t - c, -h1, h1);
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
}

bool boxBox(const Shape& sa, const Transform& ta, const Shape& sb, const Transform& tb, float margin,
            ContactManifold& m) {
    BoxFrame A = frameOf(sa, ta), B = frameOf(sb, tb);
    Vec3 t = B.c - A.c;

    float sepA = -1e30f, sepB = -1e30f, sepE = -1e30f;
    int axA = 0, axB = 0, edgeI = 0, edgeJ = 0;
    Vec3 nA, nB, nE;

    for (int i = 0; i < 3; ++i) {
        const Vec3& L = A.axis[i];
        float d = dot(t, L);
        float s = std::fabs(d) - (A.e[i] + projectRadius(B, L));
        if (s > margin) return false;
        if (s > sepA) {
            sepA = s;
            axA = i;
            nA = d >= 0 ? L : -L;
        }
    }
    for (int j = 0; j < 3; ++j) {
        const Vec3& L = B.axis[j];
        float d = dot(t, L);
        float s = std::fabs(d) - (projectRadius(A, L) + B.e[j]);
        if (s > margin) return false;
        if (s > sepB) {
            sepB = s;
            axB = j;
            nB = d >= 0 ? L : -L;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Vec3 L = cross(A.axis[i], B.axis[j]);
            float len = length(L);
            if (len < 1e-4f) continue;  // parallel edges: covered by the face axes
            L = L / len;
            float d = dot(t, L);
            float s = std::fabs(d) - (projectRadius(A, L) + projectRadius(B, L));
            if (s > margin) return false;
            if (s > sepE) {
                sepE = s;
                edgeI = i;
                edgeJ = j;
                nE = d >= 0 ? L : -L;
            }
        }
    }

    // Prefer faces (A, then B) unless another axis is clearly better; stops the manifold flickering
    // between nearly equal axes from frame to frame.
    auto clearlyBetter = [](float x, float ref) { return x > ref + 0.05f * std::fabs(ref) + 0.005f; };
    bool useB = clearlyBetter(sepB, sepA);
    float faceSep = useB ? sepB : sepA;

    if (clearlyBetter(sepE, faceSep)) {
        // Edge-edge: one contact point between the two support edges.
        const Vec3& dA = A.axis[edgeI];
        const Vec3& dB = B.axis[edgeJ];
        Vec3 cA = A.c, cB = B.c;
        for (int k = 0; k < 3; ++k) {
            if (k != edgeI) cA += A.axis[k] * (A.e[k] * (dot(A.axis[k], nE) > 0 ? 1.0f : -1.0f));
            if (k != edgeJ) cB += B.axis[k] * (B.e[k] * (dot(B.axis[k], nE) > 0 ? -1.0f : 1.0f));
        }
        Vec3 pA, pB;
        closestPointsSegments(cA, dA, A.e[edgeI], cB, dB, B.e[edgeJ], pA, pB);
        m.normal = nE;
        m.count = 1;
        m.points[0] = {(pA + pB) * 0.5f, -sepE, 0x20000u | static_cast<uint32_t>(edgeI * 3 + edgeJ)};
        return true;
    }

    if (useB)
        faceContact(B, axB, -nB, A, margin, 0, m, true);
    else
        faceContact(A, axA, nA, B, margin, 0, m, false);
    return m.count > 0;
}

}  // namespace

bool collide(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb, float margin,
             ContactManifold& out) {
    out.count = 0;
    if (a.type == ShapeType::Sphere && b.type == ShapeType::Sphere) return sphereSphere(a, ta, b, tb, margin, out);
    if (a.type == ShapeType::Box && b.type == ShapeType::Box) return boxBox(a, ta, b, tb, margin, out);
    if (a.type == ShapeType::Box) return boxSphere(a, ta, b, tb, margin, out);
    // Sphere vs box: solve as box vs sphere and flip the normal.
    bool hit = boxSphere(b, tb, a, ta, margin, out);
    out.normal = -out.normal;
    return hit;
}

}  // namespace rbe
