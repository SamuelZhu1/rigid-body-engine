#include <gtest/gtest.h>

#include <cmath>

#include "rbe/world.hpp"
#include "scenes.hpp"

using namespace rbe;

namespace {
constexpr float kDt = 1.0f / 60.0f;

WorldSettings noSleep() {
    WorldSettings s;
    s.allowSleep = false;
    return s;
}
}  // namespace

TEST(Physics, ProjectileFollowsSemiImplicitEuler) {
    WorldSettings s;
    World w(s);
    BodyDef d;
    d.shape.type = ShapeType::Sphere;
    d.position = {0, 100, 0};
    d.linearVelocity = {3, 10, 0};
    BodyId b = w.createBody(d);
    const int n = 120;
    for (int i = 0; i < n; ++i) w.step(kDt);
    // Closed form of v += g h; x += v h after n steps.
    float g = 9.81f, h = kDt;
    float y = 100 + 10 * n * h - g * h * h * n * (n + 1) / 2.0f;
    EXPECT_NEAR(w.position(b).x, 3 * n * h, 1e-4f);
    EXPECT_NEAR(w.position(b).y, y, 2e-3f);  // float accumulation over 120 steps at y ~ 100
}

TEST(Physics, RestingBoxStaysPut) {
    World w(noSleep());
    scenes::addGround(w);
    BodyDef d;
    d.position = {0, 0.5f, 0};
    BodyId b = w.createBody(d);
    for (int i = 0; i < 600; ++i) w.step(kDt);
    EXPECT_LT(length(w.position(b) - Vec3{0, 0.5f, 0}), 0.001f + 0.005f);  // within slop + 1 mm
    EXPECT_LT(std::fabs(w.position(b).x), 0.001f);
    EXPECT_LT(std::fabs(w.position(b).z), 0.001f);
}

TEST(Physics, TwentyBoxStackStandsForTenMinutes) {
    World w(noSleep());
    scenes::addGround(w);
    auto ids = scenes::stack(w, 20);
    // Let it settle for 2 s, then the top box must move less than 5 cm over 10 simulated minutes.
    for (int i = 0; i < 120; ++i) w.step(kDt);
    Vec3 start = w.position(ids.back());
    for (int i = 0; i < 36000; ++i) w.step(kDt);
    EXPECT_LT(length(w.position(ids.back()) - start), 0.05f);
}

TEST(Physics, InelasticHeadOnCollisionConservesMomentum) {
    WorldSettings s;
    s.gravity = {0, 0, 0};
    World w(s);
    BodyDef d;
    d.shape.type = ShapeType::Sphere;
    d.friction = 0;
    d.position = {-2, 0, 0};
    d.linearVelocity = {4, 0, 0};
    d.mass = 2;
    BodyId a = w.createBody(d);
    d.position = {2, 0, 0};
    d.linearVelocity = {-1, 0, 0};
    d.mass = 1;
    BodyId b = w.createBody(d);
    float p0 = 2 * 4 + 1 * -1;
    for (int i = 0; i < 120; ++i) w.step(kDt);
    float p1 = 2 * w.linearVelocity(a).x + 1 * w.linearVelocity(b).x;
    EXPECT_NEAR(p1, p0, 0.001f * std::fabs(p0));
    // Restitution 0: they end up moving together.
    EXPECT_NEAR(w.linearVelocity(a).x, w.linearVelocity(b).x, 0.01f);
}

TEST(Physics, ImplicitGyroNeverGainsEnergy) {
    WorldSettings s;
    s.gyro = GyroMode::Implicit;
    World w(s);
    scenes::spinningBox(w);
    double e0 = w.kineticEnergy(), prev = e0;
    for (int i = 0; i < 3600; ++i) {
        w.step(kDt);
        double e = w.kineticEnergy();
        ASSERT_LE(e, prev * (1 + 1e-5)) << "step " << i;
        prev = e;
    }
    EXPECT_GT(prev, 0.5 * e0);  // loses some energy, but not most of it
}

TEST(Physics, ExplicitGyroGainsEnergy) {
    WorldSettings s;
    s.gyro = GyroMode::Explicit;
    World w(s);
    scenes::spinningBox(w);
    double e0 = w.kineticEnergy();
    for (int i = 0; i < 3600; ++i) w.step(kDt);
    EXPECT_GT(w.kineticEnergy(), 1.5 * e0);
}

TEST(Physics, SameInputsGiveIdenticalState) {
    auto run = [] {
        World w;
        scenes::build(w, "pile", 400);
        for (int i = 0; i < 300; ++i) w.step(kDt);
        return w.stateHash();
    };
    uint64_t a = run(), b = run();
    EXPECT_EQ(a, b);
}

TEST(Physics, PyramidSleepsAndWakesWhenHit) {
    World w;
    scenes::addGround(w);
    auto ids = scenes::pyramid(w, 6);
    for (int i = 0; i < 600; ++i) w.step(kDt);
    for (BodyId id : ids) EXPECT_TRUE(w.isSleeping(id));
    // Drop a sphere on top.
    BodyDef d;
    d.shape.type = ShapeType::Sphere;
    d.position = {0, 7.5f, 0};
    d.linearVelocity = {0, -5, 0};
    w.createBody(d);
    for (int i = 0; i < 30; ++i) w.step(kDt);
    EXPECT_FALSE(w.isSleeping(ids.back()));
}

TEST(Physics, DominoesAllFall) {
    World w;
    scenes::addGround(w);
    auto ids = scenes::dominoes(w, 15);
    for (int i = 0; i < 600; ++i) w.step(kDt);
    // A fallen domino's long (y) axis is no longer vertical.
    Vec3 up = rotate(w.orientation(ids.back()), {0, 1, 0});
    EXPECT_LT(up.y, 0.7f);
}

TEST(Physics, SpheresComeToRest) {
    World w;
    scenes::addGround(w);
    auto ids = scenes::sphereRain(w, 50);
    for (int i = 0; i < 1200; ++i) w.step(kDt);
    for (BodyId id : ids) {
        EXPECT_GT(w.position(id).y, w.shape(id).radius - 0.02f);
        EXPECT_LT(length(w.linearVelocity(id)), 0.5f);
    }
}
