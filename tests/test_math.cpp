#include <gtest/gtest.h>

#include "rbe/math.hpp"

using namespace rbe;

namespace {
void expectNear(const Vec3& a, const Vec3& b, float tol = 1e-5f) {
    EXPECT_NEAR(a.x, b.x, tol);
    EXPECT_NEAR(a.y, b.y, tol);
    EXPECT_NEAR(a.z, b.z, tol);
}
}  // namespace

TEST(Math, CrossProductIsRightHanded) {
    expectNear(cross({1, 0, 0}, {0, 1, 0}), {0, 0, 1});
    expectNear(cross({0, 1, 0}, {0, 0, 1}), {1, 0, 0});
}

TEST(Math, SkewMatrixMatchesCross) {
    Vec3 a{0.3f, -1.2f, 2.0f}, b{-0.7f, 0.4f, 1.1f};
    expectNear(Mat3::skew(a) * b, cross(a, b));
}

TEST(Math, InverseTimesMatrixIsIdentity) {
    Mat3 m;
    m.r[0] = {2, 1, 0};
    m.r[1] = {0.5f, 3, 1};
    m.r[2] = {0, -1, 4};
    Mat3 p = inverse(m) * m;
    for (int i = 0; i < 3; ++i) expectNear(p.r[i], Mat3::identity().r[i]);
}

TEST(Math, QuaternionRotationMatchesMatrix) {
    Quat q = Quat::fromAxisAngle({1, 2, 3}, 0.9f);
    Vec3 v{0.4f, -2.0f, 1.5f};
    expectNear(rotate(q, v), toMat3(q) * v);
    expectNear(rotateInv(q, rotate(q, v)), v);
}

TEST(Math, QuarterTurnAboutZ) {
    Quat q = Quat::fromAxisAngle({0, 0, 1}, 3.14159265f / 2);
    expectNear(rotate(q, {1, 0, 0}), {0, 1, 0});
}

TEST(Math, QuaternionProductComposesRotations) {
    Quat a = Quat::fromAxisAngle({0, 1, 0}, 0.5f), b = Quat::fromAxisAngle({1, 0, 0}, -1.1f);
    Vec3 v{1, 2, 3};
    expectNear(rotate(a * b, v), rotate(a, rotate(b, v)));
}

TEST(Math, IntegratingSpinMatchesAxisAngle) {
    // 1000 small steps of 1 rad/s about y for 1 s should land near a 1 rad rotation.
    Quat q;
    for (int i = 0; i < 1000; ++i) q = integrate(q, {0, 1, 0}, 0.001f);
    Quat ref = Quat::fromAxisAngle({0, 1, 0}, 1.0f);
    expectNear(rotate(q, {1, 0, 0}), rotate(ref, {1, 0, 0}), 1e-3f);
}

TEST(Math, TangentBasisIsOrthonormal) {
    for (Vec3 n : {Vec3{0, 1, 0}, Vec3{1, 0, 0}, normalize(Vec3{1, 1, 1}), normalize(Vec3{-0.2f, 0.1f, 0.97f})}) {
        Vec3 t1, t2;
        tangentBasis(n, t1, t2);
        EXPECT_NEAR(dot(n, t1), 0, 1e-5f);
        EXPECT_NEAR(dot(n, t2), 0, 1e-5f);
        EXPECT_NEAR(dot(t1, t2), 0, 1e-5f);
        EXPECT_NEAR(length(t1), 1, 1e-5f);
        EXPECT_NEAR(length(t2), 1, 1e-5f);
    }
}

TEST(Math, AabbOverlapAndMerge) {
    AABB a{{0, 0, 0}, {1, 1, 1}}, b{{0.5f, 0.5f, 0.5f}, {2, 2, 2}}, c{{1.5f, 0, 0}, {2, 1, 1}};
    EXPECT_TRUE(a.overlaps(b));
    EXPECT_FALSE(a.overlaps(c));
    AABB m = merge(a, c);
    EXPECT_TRUE(m.contains(a));
    EXPECT_TRUE(m.contains(c));
    EXPECT_FLOAT_EQ(a.surfaceArea(), 6.0f);
}
