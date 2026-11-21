#include <gtest/gtest.h>

#include "collision.hpp"

using namespace rbe;

namespace {
Transform tf(const Vec3& p, const Quat& q = Quat::identity()) { return {p, q, toMat3(q)}; }
Shape box(const Vec3& h) {
    Shape s;
    s.type = ShapeType::Box;
    s.halfExtents = h;
    return s;
}
Shape sphere(float r) {
    Shape s;
    s.type = ShapeType::Sphere;
    s.radius = r;
    return s;
}
}  // namespace

TEST(Collision, SpheresOverlapping) {
    ContactManifold m;
    ASSERT_TRUE(collide(sphere(1), tf({0, 0, 0}), sphere(1), tf({1.5f, 0, 0}), 0, m));
    EXPECT_EQ(m.count, 1);
    EXPECT_NEAR(m.normal.x, 1, 1e-6f);
    EXPECT_NEAR(m.points[0].depth, 0.5f, 1e-6f);
    EXPECT_NEAR(m.points[0].position.x, 0.75f, 1e-6f);
}

TEST(Collision, SpheresSeparated) {
    ContactManifold m;
    EXPECT_FALSE(collide(sphere(1), tf({0, 0, 0}), sphere(1), tf({2.5f, 0, 0}), 0, m));
}

TEST(Collision, SphereRestingOnBoxFace) {
    ContactManifold m;
    ASSERT_TRUE(collide(box({1, 1, 1}), tf({0, 0, 0}), sphere(0.5f), tf({0.2f, 1.4f, 0}), 0, m));
    EXPECT_NEAR(m.normal.y, 1, 1e-6f);
    EXPECT_NEAR(m.points[0].depth, 0.1f, 1e-5f);
}

TEST(Collision, SphereOnBoxEdgeHasDiagonalNormal) {
    ContactManifold m;
    ASSERT_TRUE(collide(box({1, 1, 1}), tf({0, 0, 0}), sphere(0.5f), tf({1.3f, 1.3f, 0}), 0, m));
    EXPECT_NEAR(m.normal.x, m.normal.y, 1e-5f);
    EXPECT_NEAR(m.normal.x, 0.70710678f, 1e-5f);
}

TEST(Collision, SphereVsBoxFlipsNormal) {
    ContactManifold m;
    ASSERT_TRUE(collide(sphere(0.5f), tf({0, 1.4f, 0}), box({1, 1, 1}), tf({0, 0, 0}), 0, m));
    EXPECT_NEAR(m.normal.y, -1, 1e-6f);  // from the sphere toward the box
}

TEST(Collision, SphereCentreInsideBox) {
    ContactManifold m;
    ASSERT_TRUE(collide(box({1, 1, 1}), tf({0, 0, 0}), sphere(0.2f), tf({0, 0.9f, 0}), 0, m));
    EXPECT_NEAR(m.normal.y, 1, 1e-6f);
    EXPECT_NEAR(m.points[0].depth, 0.3f, 1e-5f);
}

TEST(Collision, BoxRestingFaceToFaceGivesFourPoints) {
    ContactManifold m;
    ASSERT_TRUE(collide(box({5, 0.5f, 5}), tf({0, -0.5f, 0}), box({0.5f, 0.5f, 0.5f}), tf({0, 0.49f, 0}), 0, m));
    EXPECT_EQ(m.count, 4);
    EXPECT_NEAR(m.normal.y, 1, 1e-6f);
    for (int i = 0; i < m.count; ++i) {
        EXPECT_NEAR(m.points[i].depth, 0.01f, 1e-5f);
        EXPECT_NEAR(std::fabs(m.points[i].position.x), 0.5f, 1e-5f);
        EXPECT_NEAR(std::fabs(m.points[i].position.z), 0.5f, 1e-5f);
    }
}

TEST(Collision, NormalPointsFromAToBWhenBIsReference) {
    // Small box A sitting under a big box B: the reference face comes from B, normal must still point A -> B.
    ContactManifold m;
    ASSERT_TRUE(collide(box({0.5f, 0.5f, 0.5f}), tf({0, 0.49f, 0}), box({5, 0.5f, 5}), tf({0, 1.48f, 0}), 0, m));
    EXPECT_NEAR(m.normal.y, 1, 1e-6f);
    EXPECT_EQ(m.count, 4);
}

TEST(Collision, RotatedBoxOnGroundGivesFourPointsInsideFace) {
    ContactManifold m;
    Quat q = Quat::fromAxisAngle({0, 1, 0}, 0.785398f);  // 45 deg about the vertical
    ASSERT_TRUE(collide(box({5, 0.5f, 5}), tf({0, -0.5f, 0}), box({0.5f, 0.5f, 0.5f}), tf({0, 0.48f, 0}, q), 0, m));
    EXPECT_EQ(m.count, 4);
    for (int i = 0; i < m.count; ++i) EXPECT_NEAR(m.points[i].depth, 0.02f, 1e-5f);
}

TEST(Collision, SmallBoxOnBigBoxClipsToSmallFace) {
    // Big box resting on a small one: incident face is bigger and must be clipped to the reference face.
    ContactManifold m;
    ASSERT_TRUE(collide(box({0.5f, 0.5f, 0.5f}), tf({0, 0, 0}), box({2, 0.5f, 2}), tf({0, 0.99f, 0}), 0, m));
    EXPECT_EQ(m.count, 4);
    for (int i = 0; i < m.count; ++i) {
        EXPECT_LE(std::fabs(m.points[i].position.x), 0.5f + 1e-5f);
        EXPECT_LE(std::fabs(m.points[i].position.z), 0.5f + 1e-5f);
    }
}

TEST(Collision, CrossedEdgesGiveOnePoint) {
    // Two long thin boxes crossing at right angles, each turned 45 deg about its long axis so edges meet.
    Quat qa = Quat::fromAxisAngle({1, 0, 0}, 0.785398f);
    Quat qb = Quat::fromAxisAngle({0, 1, 0}, 1.570796f) * Quat::fromAxisAngle({1, 0, 0}, 0.785398f);
    float d = 2 * 0.70710678f - 0.01f;  // each diamond tip sticks out 0.5 * sqrt(2) from its centre
    ContactManifold m;
    ASSERT_TRUE(collide(box({2, 0.5f, 0.5f}), tf({0, 0, 0}, qa), box({2, 0.5f, 0.5f}), tf({0, d, 0}, qb), 0, m));
    EXPECT_EQ(m.count, 1);
    EXPECT_NEAR(m.normal.y, 1, 1e-4f);
    EXPECT_NEAR(m.points[0].depth, 0.01f, 1e-4f);
    EXPECT_NEAR(m.points[0].position.x, 0, 1e-4f);
    EXPECT_NEAR(m.points[0].position.z, 0, 1e-4f);
}

TEST(Collision, SeparatedBoxesReportNothing) {
    ContactManifold m;
    EXPECT_FALSE(collide(box({0.5f, 0.5f, 0.5f}), tf({0, 0, 0}), box({0.5f, 0.5f, 0.5f}), tf({1.2f, 0.3f, 0}), 0, m));
}

TEST(Collision, FeatureIdsStableUnderSmallMotion) {
    // Warm starting matches points by id, so a tiny slide must not change the ids.
    ContactManifold m1, m2;
    Shape g = box({5, 0.5f, 5}), b = box({0.5f, 0.5f, 0.5f});
    ASSERT_TRUE(collide(g, tf({0, -0.5f, 0}), b, tf({0, 0.49f, 0}), 0, m1));
    ASSERT_TRUE(collide(g, tf({0, -0.5f, 0}), b, tf({0.003f, 0.491f, -0.002f}), 0, m2));
    ASSERT_EQ(m1.count, m2.count);
    for (int i = 0; i < m1.count; ++i) {
        bool found = false;
        for (int j = 0; j < m2.count; ++j) found |= m1.points[i].id == m2.points[j].id;
        EXPECT_TRUE(found);
    }
}

TEST(Collision, BoxAabbCoversRotatedBox) {
    Quat q = Quat::fromAxisAngle({0, 0, 1}, 0.785398f);
    AABB a = computeAABB(box({1, 1, 1}), tf({0, 0, 0}, q));
    EXPECT_NEAR(a.hi.x, 1.41421356f, 1e-5f);
    EXPECT_NEAR(a.hi.z, 1.0f, 1e-5f);
}
