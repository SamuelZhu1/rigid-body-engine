// Public API of the rigid body engine.
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "rbe/math.hpp"

namespace rbe {

enum class ShapeType : uint8_t { Sphere, Box };

struct Shape {
    ShapeType type = ShapeType::Box;
    Vec3 halfExtents{0.5f, 0.5f, 0.5f};  // box
    float radius = 0.5f;                 // sphere
};

enum class GyroMode : uint8_t {
    None,      // ignore the gyroscopic term (many game engines do)
    Explicit,  // w += h * I^-1 * -(w x Iw); gains energy, blows up fast spinners
    Implicit,  // one Newton step of backward Euler (Catto, GDC 2015)
};

enum class BroadphaseMode : uint8_t { Tree, BruteForce };

struct WorldSettings {
    Vec3 gravity{0.0f, -9.81f, 0.0f};
    // Each step is split into substeps; every substep integrates, warm starts, solves contacts with
    // soft constraints and then relaxes (Box2D v3 style "soft step").
    int substeps = 4;
    float contactHertz = 30.0f;        // stiffness of the soft contact spring that pushes overlap out
    float contactDampingRatio = 10.0f;  // heavily damped so overlap is removed without bounce
    bool warmStarting = true;
    bool allowSleep = true;
    GyroMode gyro = GyroMode::Implicit;
    BroadphaseMode broadphase = BroadphaseMode::Tree;
    float linearSlop = 0.005f;        // penetration allowed before position correction kicks in
    float maxCorrectionSpeed = 3.0f;  // m/s cap on push-out speed so deep overlaps don't explode
    float aabbMargin = 0.1f;          // fat AABB margin in the broadphase tree
    float sleepTime = 0.5f;
    float sleepLinear = 0.05f;
    float sleepAngular = 2.0f * 3.14159265f / 180.0f;
};

struct BodyDef {
    Shape shape;
    Vec3 position{0, 0, 0};
    Quat orientation = Quat::identity();
    Vec3 linearVelocity{0, 0, 0};
    Vec3 angularVelocity{0, 0, 0};
    float mass = 1.0f;  // 0 = static
    float friction = 0.6f;
};

struct BodyId {
    uint32_t index = UINT32_MAX;
    bool valid() const { return index != UINT32_MAX; }
};

struct StepStats {
    int pairs = 0;          // broadphase candidate pairs
    int manifolds = 0;      // touching pairs
    int contactPoints = 0;  // total contact points
    int awakeBodies = 0;
    double msIntegrate = 0, msBroadphase = 0, msNarrowphase = 0, msSolve = 0, msIslands = 0, msTotal = 0;
};

// Read-only view of one contact point, for debug drawing.
struct DebugContact {
    Vec3 point, normal;
    float depth;
};

class World {
public:
    explicit World(const WorldSettings& s = {});
    ~World();
    World(const World&) = delete;
    World& operator=(const World&) = delete;

    BodyId createBody(const BodyDef& def);
    void step(float dt);

    // Body state.
    int bodyCount() const;
    Vec3 position(BodyId id) const;
    Quat orientation(BodyId id) const;
    Vec3 linearVelocity(BodyId id) const;
    Vec3 angularVelocity(BodyId id) const;
    const Shape& shape(BodyId id) const;
    bool isStatic(BodyId id) const;
    bool isSleeping(BodyId id) const;
    void setVelocity(BodyId id, const Vec3& v, const Vec3& w);
    void applyImpulse(BodyId id, const Vec3& impulse, const Vec3& worldPoint);

    // Kinetic energy (linear + rotational) of all dynamic bodies.
    double kineticEnergy() const;
    // Angular momentum of one body about its centre, world frame.
    Vec3 angularMomentum(BodyId id) const;
    // FNV-1a over the raw bytes of every position and orientation. Same build + same inputs => same hash.
    uint64_t stateHash() const;

    WorldSettings& settings();
    const StepStats& stats() const;
    std::vector<DebugContact> debugContacts() const;
    std::vector<AABB> debugTreeNodes() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace rbe
