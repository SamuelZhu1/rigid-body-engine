// Small 3D math library: vectors, quaternions, 3x3 matrices, AABBs.
// Written from scratch on purpose; float everywhere, no SIMD.
#pragma once

#include <cmath>
#include <cstdint>

namespace rbe {

struct Vec3 {
    float x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    float& operator[](int i) { return (&x)[i]; }
    float operator[](int i) const { return (&x)[i]; }

    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
    Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3& a, float s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(float s, const Vec3& a) { return a * s; }
inline Vec3 operator/(const Vec3& a, float s) { return a * (1.0f / s); }

inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float lengthSq(const Vec3& a) { return dot(a, a); }
inline float length(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(const Vec3& a) {
    float l = length(a);
    return l > 1e-12f ? a / l : Vec3{0, 0, 0};
}
inline Vec3 vmin(const Vec3& a, const Vec3& b) {
    return {std::fmin(a.x, b.x), std::fmin(a.y, b.y), std::fmin(a.z, b.z)};
}
inline Vec3 vmax(const Vec3& a, const Vec3& b) {
    return {std::fmax(a.x, b.x), std::fmax(a.y, b.y), std::fmax(a.z, b.z)};
}
inline Vec3 vabs(const Vec3& a) { return {std::fabs(a.x), std::fabs(a.y), std::fabs(a.z)}; }

// Builds two unit vectors perpendicular to n (and each other). n must be unit length.
inline void tangentBasis(const Vec3& n, Vec3& t1, Vec3& t2) {
    // Pick the axis least aligned with n to avoid a degenerate cross product.
    if (std::fabs(n.x) >= 0.57735f)
        t1 = normalize(Vec3{n.y, -n.x, 0.0f});
    else
        t1 = normalize(Vec3{0.0f, n.z, -n.y});
    t2 = cross(n, t1);
}

// Row-major 3x3 matrix.
struct Mat3 {
    Vec3 r[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    static Mat3 identity() { return {}; }
    static Mat3 zero() {
        Mat3 m;
        m.r[0] = m.r[1] = m.r[2] = {0, 0, 0};
        return m;
    }
    static Mat3 diagonal(const Vec3& d) {
        Mat3 m = zero();
        m.r[0].x = d.x;
        m.r[1].y = d.y;
        m.r[2].z = d.z;
        return m;
    }
    // Matrix form of the cross product: skew(a) * b == cross(a, b).
    static Mat3 skew(const Vec3& a) {
        Mat3 m;
        m.r[0] = {0, -a.z, a.y};
        m.r[1] = {a.z, 0, -a.x};
        m.r[2] = {-a.y, a.x, 0};
        return m;
    }

    Vec3 col(int i) const { return {r[0][i], r[1][i], r[2][i]}; }
};

inline Vec3 operator*(const Mat3& m, const Vec3& v) { return {dot(m.r[0], v), dot(m.r[1], v), dot(m.r[2], v)}; }
inline Mat3 transpose(const Mat3& m) {
    Mat3 t;
    t.r[0] = m.col(0);
    t.r[1] = m.col(1);
    t.r[2] = m.col(2);
    return t;
}
inline Mat3 operator*(const Mat3& a, const Mat3& b) {
    Mat3 out;
    Mat3 bt = transpose(b);
    for (int i = 0; i < 3; ++i) out.r[i] = {dot(a.r[i], bt.r[0]), dot(a.r[i], bt.r[1]), dot(a.r[i], bt.r[2])};
    return out;
}
inline Mat3 operator+(const Mat3& a, const Mat3& b) {
    Mat3 out;
    for (int i = 0; i < 3; ++i) out.r[i] = a.r[i] + b.r[i];
    return out;
}
inline Mat3 operator-(const Mat3& a, const Mat3& b) {
    Mat3 out;
    for (int i = 0; i < 3; ++i) out.r[i] = a.r[i] - b.r[i];
    return out;
}
inline Mat3 operator*(const Mat3& a, float s) {
    Mat3 out;
    for (int i = 0; i < 3; ++i) out.r[i] = a.r[i] * s;
    return out;
}
inline float determinant(const Mat3& m) { return dot(m.r[0], cross(m.r[1], m.r[2])); }
inline Mat3 inverse(const Mat3& m) {
    // Rows of the inverse transpose are the cross products of the rows.
    float det = determinant(m);
    float inv = std::fabs(det) > 1e-20f ? 1.0f / det : 0.0f;
    Mat3 adjT;
    adjT.r[0] = cross(m.r[1], m.r[2]) * inv;
    adjT.r[1] = cross(m.r[2], m.r[0]) * inv;
    adjT.r[2] = cross(m.r[0], m.r[1]) * inv;
    return transpose(adjT);
}

// Unit quaternion rotation, w + xi + yj + zk.
struct Quat {
    float w = 1, x = 0, y = 0, z = 0;

    static Quat identity() { return {}; }
    static Quat fromAxisAngle(const Vec3& axis, float angle) {
        Vec3 a = normalize(axis);
        float s = std::sin(0.5f * angle);
        return {std::cos(0.5f * angle), a.x * s, a.y * s, a.z * s};
    }
    Vec3 v() const { return {x, y, z}; }
};

inline Quat operator*(const Quat& a, const Quat& b) {
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}
inline Quat conjugate(const Quat& q) { return {q.w, -q.x, -q.y, -q.z}; }
inline Quat normalize(const Quat& q) {
    float l = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    float inv = l > 1e-12f ? 1.0f / l : 0.0f;
    return {q.w * inv, q.x * inv, q.y * inv, q.z * inv};
}
inline Vec3 rotate(const Quat& q, const Vec3& v) {
    // v' = v + 2w(u x v) + 2 u x (u x v), u = vector part.
    Vec3 u = q.v();
    Vec3 t = cross(u, v) * 2.0f;
    return v + t * q.w + cross(u, t);
}
inline Vec3 rotateInv(const Quat& q, const Vec3& v) { return rotate(conjugate(q), v); }
inline Mat3 toMat3(const Quat& q) {
    float xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    float wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    Mat3 m;
    m.r[0] = {1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy)};
    m.r[1] = {2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx)};
    m.r[2] = {2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)};
    return m;
}
// First-order update of an orientation by world angular velocity w over h seconds:
// q' = q + h/2 * (0, w) * q, then renormalize.
inline Quat integrate(const Quat& q, const Vec3& w, float h) {
    Quat spin{0, w.x, w.y, w.z};
    Quat d = spin * q;
    float k = 0.5f * h;
    return normalize(Quat{q.w + d.w * k, q.x + d.x * k, q.y + d.y * k, q.z + d.z * k});
}

struct AABB {
    Vec3 lo{0, 0, 0}, hi{0, 0, 0};

    bool overlaps(const AABB& b) const {
        return lo.x <= b.hi.x && hi.x >= b.lo.x && lo.y <= b.hi.y && hi.y >= b.lo.y && lo.z <= b.hi.z &&
               hi.z >= b.lo.z;
    }
    bool contains(const AABB& b) const {
        return lo.x <= b.lo.x && lo.y <= b.lo.y && lo.z <= b.lo.z && hi.x >= b.hi.x && hi.y >= b.hi.y &&
               hi.z >= b.hi.z;
    }
    float surfaceArea() const {
        Vec3 d = hi - lo;
        return 2.0f * (d.x * d.y + d.y * d.z + d.z * d.x);
    }
    AABB fattened(float m) const { return {lo - Vec3{m, m, m}, hi + Vec3{m, m, m}}; }
};

inline AABB merge(const AABB& a, const AABB& b) { return {vmin(a.lo, b.lo), vmax(a.hi, b.hi)}; }

}  // namespace rbe
