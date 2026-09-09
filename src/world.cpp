#include "rbe/world.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <numeric>
#include <unordered_map>
#include <utility>

#include "collision.hpp"
#include "dynamic_tree.hpp"
#include "profile.hpp"

namespace rbe {

namespace {

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

constexpr float kContactMargin = 0.02f;  // keep touching contacts alive slightly before overlap (speculative)

struct SolverPoint {
    Vec3 rA, rB;               // contact point relative to each body at the start of the step
    Vec3 localA, localB;       // the same anchors in each body's local frame
    float baseSeparation = 0;  // -penetration at the start of the step
    float normalImpulse = 0, tangentImpulse[2] = {0, 0};
    float normalMass = 0, tangentMass[2] = {0, 0};
    uint32_t id = 0;
    Vec3 position;
};

struct SolverManifold {
    uint32_t a = 0, b = 0;
    Vec3 normal, t1, t2;
    float friction = 0.6f;
    int count = 0;
    SolverPoint pts[8];
};

// Soft constraint coefficients (Box2D v3 "soft step"): a damped spring of the given stiffness pushes
// overlap out over several substeps instead of all at once.
struct Softness {
    float biasRate = 0, massScale = 1, impulseScale = 0;
};
Softness makeSoft(float hertz, float zeta, float h) {
    if (hertz == 0.0f) return {};
    float omega = 2.0f * 3.14159265f * hertz;
    float a1 = 2.0f * zeta + h * omega;
    float a2 = h * omega * a1;
    float a3 = 1.0f / (1.0f + a2);
    return {omega / a1, a2 * a3, a3};
}

uint64_t pairKey(uint32_t a, uint32_t b) { return (static_cast<uint64_t>(a) << 32) | b; }

Vec3 boxInertia(float m, const Vec3& e) {
    return {m / 3.0f * (e.y * e.y + e.z * e.z), m / 3.0f * (e.x * e.x + e.z * e.z), m / 3.0f * (e.x * e.x + e.y * e.y)};
}

}  // namespace

struct World::Impl {
    WorldSettings s;

    // Bodies, structure of arrays: the hot loops only touch the arrays they need.
    std::vector<Vec3> pos, vel, angVel;
    std::vector<Quat> rot;
    std::vector<Mat3> rotM, invInertiaWorld;
    std::vector<Vec3> inertiaLocal, invInertiaLocal;
    std::vector<float> invMass, mass, friction, sleepTimer;
    std::vector<uint8_t> isStatic, sleeping;
    std::vector<Shape> shapes;
    std::vector<int32_t> proxy;

    DynamicTree tree;
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    std::vector<SolverManifold> manifolds;
    std::unordered_map<uint64_t, SolverManifold> cache;  // last step's manifolds, for warm starting
    std::vector<uint32_t> parent;                        // union-find for islands
    StepStats stats;

    Transform transform(uint32_t i) const { return {pos[i], rot[i], rotM[i]}; }
    bool awakeDynamic(uint32_t i) const { return !isStatic[i] && !sleeping[i]; }

    void updateInertia(uint32_t i) {
        rotM[i] = toMat3(rot[i]);
        invInertiaWorld[i] = rotM[i] * Mat3::diagonal(invInertiaLocal[i]) * transpose(rotM[i]);
    }

    void integrateVelocities(float h);
    void broadphase();
    void narrowphase(float h);
    void warmStart();
    void solveContacts(float h, bool useBias);
    void integratePositions(float h);
    void islandsAndSleep(float h);

    void applyImpulse(uint32_t a, uint32_t b, const Vec3& rA, const Vec3& rB, const Vec3& P) {
        vel[a] -= P * invMass[a];
        angVel[a] -= invInertiaWorld[a] * cross(rA, P);
        vel[b] += P * invMass[b];
        angVel[b] += invInertiaWorld[b] * cross(rB, P);
    }

    uint32_t find(uint32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
};

World::World(const WorldSettings& s) : impl_(std::make_unique<Impl>()) { impl_->s = s; }
World::~World() = default;

BodyId World::createBody(const BodyDef& d) {
    Impl& w = *impl_;
    uint32_t i = static_cast<uint32_t>(w.pos.size());
    bool stat = d.mass <= 0.0f;
    w.pos.push_back(d.position);
    w.rot.push_back(normalize(d.orientation));
    w.vel.push_back(stat ? Vec3{} : d.linearVelocity);
    w.angVel.push_back(stat ? Vec3{} : d.angularVelocity);
    w.rotM.emplace_back();
    w.invInertiaWorld.emplace_back();
    w.shapes.push_back(d.shape);
    w.friction.push_back(d.friction);
    w.isStatic.push_back(stat ? 1 : 0);
    w.sleeping.push_back(0);
    w.sleepTimer.push_back(0);
    w.mass.push_back(stat ? 0.0f : d.mass);
    w.invMass.push_back(stat ? 0.0f : 1.0f / d.mass);

    Vec3 I{0, 0, 0}, invI{0, 0, 0};
    if (!stat) {
        if (d.shape.type == ShapeType::Box) {
            I = boxInertia(d.mass, d.shape.halfExtents);
        } else {
            float k = 0.4f * d.mass * d.shape.radius * d.shape.radius;
            I = {k, k, k};
        }
        invI = {1.0f / I.x, 1.0f / I.y, 1.0f / I.z};
    }
    w.inertiaLocal.push_back(I);
    w.invInertiaLocal.push_back(invI);
    w.updateInertia(i);
    w.proxy.push_back(w.tree.createProxy(computeAABB(d.shape, w.transform(i)), w.s.aabbMargin, i));
    return {i};
}

void World::Impl::integrateVelocities(float h) {
    RBE_ZONE("integrate velocities");
    const uint32_t n = static_cast<uint32_t>(pos.size());
    for (uint32_t i = 0; i < n; ++i) {
        if (!awakeDynamic(i)) continue;
        vel[i] += s.gravity * h;

        Vec3& w = angVel[i];
        switch (s.gyro) {
            case GyroMode::None:
                break;
            case GyroMode::Explicit: {
                Mat3 Iw = rotM[i] * Mat3::diagonal(inertiaLocal[i]) * transpose(rotM[i]);
                w += invInertiaWorld[i] * (-cross(w, Iw * w)) * h;
                break;
            }
            case GyroMode::Implicit: {
                // Solve Ib (w' - w) + h w' x (Ib w') = 0 in the body frame with one Newton step from w.
                const Mat3& R = rotM[i];
                Vec3 wb = transpose(R) * w;
                Mat3 Ib = Mat3::diagonal(inertiaLocal[i]);
                Vec3 Iwb = Ib * wb;
                Vec3 f = cross(wb, Iwb) * h;
                Mat3 J = Ib + (Mat3::skew(wb) * Ib - Mat3::skew(Iwb)) * h;
                wb = wb - inverse(J) * f;
                w = R * wb;
                break;
            }
        }
    }
}

void World::Impl::broadphase() {
    RBE_ZONE("broadphase");
    pairs.clear();
    const uint32_t n = static_cast<uint32_t>(pos.size());

    for (uint32_t i = 0; i < n; ++i) {
        if (!awakeDynamic(i)) continue;
        tree.moveProxy(proxy[i], computeAABB(shapes[i], transform(i)), s.aabbMargin);
    }

    if (s.broadphase == BroadphaseMode::BruteForce) {
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t j = i + 1; j < n; ++j) {
                if (!awakeDynamic(i) && !awakeDynamic(j)) continue;
                if (tree.fatAABB(proxy[i]).overlaps(tree.fatAABB(proxy[j]))) pairs.emplace_back(i, j);
            }
        }
        return;
    }

    for (uint32_t i = 0; i < n; ++i) {
        if (!awakeDynamic(i)) continue;
        tree.query(tree.fatAABB(proxy[i]), [&](uint32_t j) {
            if (j == i) return;
            // An awake dynamic partner finds this pair itself, so only the lower index keeps it.
            if (awakeDynamic(j) && j < i) return;
            pairs.emplace_back(std::min(i, j), std::max(i, j));
        });
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
}

void World::Impl::narrowphase(float h) {
    RBE_ZONE("narrowphase");
    (void)h;
    manifolds.clear();
    for (auto [a, b] : pairs) {
        ContactManifold cm;
        if (!collide(shapes[a], transform(a), shapes[b], transform(b), kContactMargin, cm)) continue;

        // A sleeping body touched by an awake one wakes up.
        if (sleeping[a]) {
            sleeping[a] = 0;
            sleepTimer[a] = 0;
        }
        if (sleeping[b]) {
            sleeping[b] = 0;
            sleepTimer[b] = 0;
        }

        SolverManifold m;
        m.a = a;
        m.b = b;
        m.normal = cm.normal;
        tangentBasis(m.normal, m.t1, m.t2);
        m.friction = std::sqrt(friction[a] * friction[b]);
        m.count = cm.count;

        const SolverManifold* old = nullptr;
        if (s.warmStarting) {
            auto it = cache.find(pairKey(a, b));
            if (it != cache.end()) old = &it->second;
        }

        auto effMass = [&](const Vec3& rAp, const Vec3& rBp, const Vec3& dir) {
            Vec3 ra = cross(rAp, dir), rb = cross(rBp, dir);
            float kk = invMass[a] + invMass[b] + dot(ra, invInertiaWorld[a] * ra) + dot(rb, invInertiaWorld[b] * rb);
            return kk > 0 ? 1.0f / kk : 0.0f;
        };
        uint32_t usedOld = 0;
        for (int k = 0; k < cm.count; ++k) {
            SolverPoint& p = m.pts[k];
            p.position = cm.points[k].position;
            p.baseSeparation = -cm.points[k].depth;
            p.id = cm.points[k].id;
            p.rA = p.position - pos[a];
            p.rB = p.position - pos[b];
            p.localA = rotateInv(rot[a], p.rA);
            p.localB = rotateInv(rot[b], p.rB);
            p.normalMass = effMass(p.rA, p.rB, m.normal);
            p.tangentMass[0] = effMass(p.rA, p.rB, m.t1);
            p.tangentMass[1] = effMass(p.rA, p.rB, m.t2);

            if (old) {
                // Match last frame's point by feature id, or failing that by position. Equal-size boxes
                // in a stack touch face to face exactly, so a tiny tilt moves a corner just outside the
                // other face and the clipper hands back a new point with a new id. Without the position
                // fallback that point loses its warm start and tall stacks shake apart. Each old point is
                // handed out once, or two new points near one corner would both get its impulse.
                int match = -1;
                for (int q = 0; q < old->count; ++q)
                    if (!(usedOld & (1u << q)) && old->pts[q].id == p.id) match = q;
                if (match < 0) {
                    float best = 0.05f * 0.05f;
                    for (int q = 0; q < old->count; ++q) {
                        if (usedOld & (1u << q)) continue;
                        float d = lengthSq(old->pts[q].localA - p.localA);
                        if (d < best) {
                            best = d;
                            match = q;
                        }
                    }
                }
                if (match >= 0) {
                    usedOld |= 1u << match;
                    const SolverPoint& op = old->pts[match];
                    p.normalImpulse = op.normalImpulse;
                    // Re-project last frame's friction impulse onto this frame's tangent basis.
                    Vec3 oldT = old->t1 * op.tangentImpulse[0] + old->t2 * op.tangentImpulse[1];
                    p.tangentImpulse[0] = dot(oldT, m.t1);
                    p.tangentImpulse[1] = dot(oldT, m.t2);
                }
            }
        }
        manifolds.push_back(m);
    }
}

void World::Impl::warmStart() {
    RBE_ZONE("warm start");
    for (SolverManifold& m : manifolds) {
        for (int k = 0; k < m.count; ++k) {
            const SolverPoint& p = m.pts[k];
            applyImpulse(m.a, m.b, p.rA, p.rB,
                         m.normal * p.normalImpulse + m.t1 * p.tangentImpulse[0] + m.t2 * p.tangentImpulse[1]);
        }
    }
}

// One pass of sequential impulses over every contact. With useBias the normal constraint is a soft
// spring that pushes out overlap; the relax pass after the position update runs it rigid with no
// push-out, which removes the velocity the push-out added (so fixing overlap doesn't make things bounce).
void World::Impl::solveContacts(float h, bool useBias) {
    RBE_ZONE("solve contacts");
    const Softness soft = makeSoft(s.contactHertz, s.contactDampingRatio, h);
    const float invH = 1.0f / h;
    for (SolverManifold& m : manifolds) {
        const uint32_t a = m.a, b = m.b;
        for (int k = 0; k < m.count; ++k) {
            SolverPoint& p = m.pts[k];
            // Current separation from where the anchors have moved to during the substeps so far.
            Vec3 d = (pos[b] + rotM[b] * p.localB) - (pos[a] + rotM[a] * p.localA);
            float sep = p.baseSeparation + dot(d, m.normal) + s.linearSlop;

            float bias = 0, massScale = 1, impulseScale = 0;
            if (sep > 0) {
                bias = sep * invH;  // speculative: allowed to close the gap this substep
            } else if (useBias) {
                bias = std::max(soft.biasRate * sep, -s.maxCorrectionSpeed);
                massScale = soft.massScale;
                impulseScale = soft.impulseScale;
            }

            Vec3 dv = vel[b] + cross(angVel[b], p.rB) - vel[a] - cross(angVel[a], p.rA);
            float vn = dot(dv, m.normal);
            float lambda = -p.normalMass * massScale * (vn + bias) - impulseScale * p.normalImpulse;
            float old = p.normalImpulse;
            p.normalImpulse = std::max(old + lambda, 0.0f);
            applyImpulse(a, b, p.rA, p.rB, m.normal * (p.normalImpulse - old));
        }
        for (int k = 0; k < m.count; ++k) {
            SolverPoint& p = m.pts[k];
            float maxF = m.friction * p.normalImpulse;
            const Vec3 tangents[2] = {m.t1, m.t2};
            for (int t = 0; t < 2; ++t) {
                Vec3 dv = vel[b] + cross(angVel[b], p.rB) - vel[a] - cross(angVel[a], p.rA);
                float lambda = -p.tangentMass[t] * dot(dv, tangents[t]);
                float old = p.tangentImpulse[t];
                p.tangentImpulse[t] = std::clamp(old + lambda, -maxF, maxF);
                applyImpulse(a, b, p.rA, p.rB, tangents[t] * (p.tangentImpulse[t] - old));
            }
        }
    }
}

void World::Impl::integratePositions(float h) {
    RBE_ZONE("integrate positions");
    const uint32_t n = static_cast<uint32_t>(pos.size());
    for (uint32_t i = 0; i < n; ++i) {
        if (!awakeDynamic(i)) continue;
        pos[i] += vel[i] * h;
        rot[i] = integrate(rot[i], angVel[i], h);
        updateInertia(i);
    }
}

void World::Impl::islandsAndSleep(float h) {
    RBE_ZONE("islands");
    const uint32_t n = static_cast<uint32_t>(pos.size());
    parent.resize(n);
    std::iota(parent.begin(), parent.end(), 0u);
    for (const SolverManifold& m : manifolds) {
        if (isStatic[m.a] || isStatic[m.b]) continue;
        uint32_t ra = find(m.a), rb = find(m.b);
        if (ra != rb) parent[std::max(ra, rb)] = std::min(ra, rb);
    }

    // Per body: how long it has been nearly still. Per island: the minimum of that.
    std::vector<float> islandMin(n, 1e30f);
    for (uint32_t i = 0; i < n; ++i) {
        if (isStatic[i]) continue;
        if (!sleeping[i]) {
            bool slow = lengthSq(vel[i]) < s.sleepLinear * s.sleepLinear &&
                        lengthSq(angVel[i]) < s.sleepAngular * s.sleepAngular;
            sleepTimer[i] = slow ? sleepTimer[i] + h : 0.0f;
        }
        uint32_t r = find(i);
        islandMin[r] = std::min(islandMin[r], sleeping[i] ? s.sleepTime : sleepTimer[i]);
    }
    int awake = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (isStatic[i]) continue;
        bool sleep = s.allowSleep && islandMin[find(i)] >= s.sleepTime;
        if (sleep && !sleeping[i]) {
            sleeping[i] = 1;
            vel[i] = {};
            angVel[i] = {};
        } else if (!sleep && sleeping[i]) {
            sleeping[i] = 0;
            sleepTimer[i] = 0;
        }
        awake += sleeping[i] ? 0 : 1;
    }
    stats.awakeBodies = awake;
}

void World::step(float h) {
    RBE_ZONE("step");
    Impl& w = *impl_;
    StepStats& st = w.stats;
    auto t0 = Clock::now();

    // Collision detection runs once per step; the solver substeps on that contact set.
    auto t = Clock::now();
    w.broadphase();
    st.msBroadphase = msSince(t);

    t = Clock::now();
    w.narrowphase(h);
    st.msNarrowphase = msSince(t);

    if (!w.s.warmStarting)
        for (SolverManifold& m : w.manifolds)
            for (int k = 0; k < m.count; ++k)
                m.pts[k].normalImpulse = m.pts[k].tangentImpulse[0] = m.pts[k].tangentImpulse[1] = 0;

    const int n = std::max(1, w.s.substeps);
    const float hs = h / static_cast<float>(n);
    double msInt = 0;
    t = Clock::now();
    for (int i = 0; i < n; ++i) {
        auto ti = Clock::now();
        w.integrateVelocities(hs);
        msInt += msSince(ti);
        w.warmStart();
        w.solveContacts(hs, true);
        w.integratePositions(hs);
        w.solveContacts(hs, false);
    }
    st.msSolve = msSince(t) - msInt;
    st.msIntegrate = msInt;

    w.cache.clear();
    int points = 0;
    for (const SolverManifold& m : w.manifolds) {
        w.cache.emplace(pairKey(m.a, m.b), m);
        points += m.count;
    }

    t = Clock::now();
    w.islandsAndSleep(h);
    st.msIslands = msSince(t);

    st.pairs = static_cast<int>(w.pairs.size());
    st.manifolds = static_cast<int>(w.manifolds.size());
    st.contactPoints = points;
    st.msTotal = msSince(t0);
}

int World::bodyCount() const { return static_cast<int>(impl_->pos.size()); }
Vec3 World::position(BodyId id) const { return impl_->pos[id.index]; }
Quat World::orientation(BodyId id) const { return impl_->rot[id.index]; }
Vec3 World::linearVelocity(BodyId id) const { return impl_->vel[id.index]; }
Vec3 World::angularVelocity(BodyId id) const { return impl_->angVel[id.index]; }
const Shape& World::shape(BodyId id) const { return impl_->shapes[id.index]; }
bool World::isStatic(BodyId id) const { return impl_->isStatic[id.index] != 0; }
bool World::isSleeping(BodyId id) const { return impl_->sleeping[id.index] != 0; }
WorldSettings& World::settings() { return impl_->s; }
const StepStats& World::stats() const { return impl_->stats; }

void World::setVelocity(BodyId id, const Vec3& v, const Vec3& w) {
    Impl& im = *impl_;
    if (im.isStatic[id.index]) return;
    im.vel[id.index] = v;
    im.angVel[id.index] = w;
    im.sleeping[id.index] = 0;
    im.sleepTimer[id.index] = 0;
}

void World::applyImpulse(BodyId id, const Vec3& P, const Vec3& point) {
    Impl& im = *impl_;
    uint32_t i = id.index;
    if (im.isStatic[i]) return;
    im.vel[i] += P * im.invMass[i];
    im.angVel[i] += im.invInertiaWorld[i] * cross(point - im.pos[i], P);
    im.sleeping[i] = 0;
    im.sleepTimer[i] = 0;
}

double World::kineticEnergy() const {
    const Impl& im = *impl_;
    double e = 0;
    for (size_t i = 0; i < im.pos.size(); ++i) {
        if (im.isStatic[i]) continue;
        Vec3 wb = transpose(im.rotM[i]) * im.angVel[i];
        const Vec3& I = im.inertiaLocal[i];
        e += 0.5 * im.mass[i] * lengthSq(im.vel[i]);
        e += 0.5 * (I.x * wb.x * wb.x + I.y * wb.y * wb.y + I.z * wb.z * wb.z);
    }
    return e;
}

Vec3 World::angularMomentum(BodyId id) const {
    const Impl& im = *impl_;
    uint32_t i = id.index;
    Vec3 wb = transpose(im.rotM[i]) * im.angVel[i];
    const Vec3& I = im.inertiaLocal[i];
    return im.rotM[i] * Vec3{I.x * wb.x, I.y * wb.y, I.z * wb.z};
}

uint64_t World::stateHash() const {
    const Impl& im = *impl_;
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](const void* data, size_t len) {
        const unsigned char* c = static_cast<const unsigned char*>(data);
        for (size_t k = 0; k < len; ++k) {
            h ^= c[k];
            h *= 1099511628211ull;
        }
    };
    mix(im.pos.data(), im.pos.size() * sizeof(Vec3));
    mix(im.rot.data(), im.rot.size() * sizeof(Quat));
    return h;
}

std::vector<DebugContact> World::debugContacts() const {
    std::vector<DebugContact> out;
    for (const SolverManifold& m : impl_->manifolds)
        for (int k = 0; k < m.count; ++k) out.push_back({m.pts[k].position, m.normal, -m.pts[k].baseSeparation});
    return out;
}

std::vector<AABB> World::debugTreeNodes() const {
    std::vector<AABB> out;
    impl_->tree.collectNodeBoxes(out);
    return out;
}

}  // namespace rbe
