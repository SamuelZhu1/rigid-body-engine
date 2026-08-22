// Narrowphase: contact manifolds for sphere and box pairs.
#pragma once

#include <cstdint>

#include "rbe/math.hpp"
#include "rbe/world.hpp"

namespace rbe {

struct Transform {
    Vec3 p;
    Quat q;
    Mat3 R;  // cached rotation matrix of q (columns are the local axes in world space)
};

struct ManifoldPoint {
    Vec3 position;     // world-space contact point (midway between the two surfaces)
    float depth = 0;   // penetration, > 0 when overlapping
    uint32_t id = 0;   // feature id, stable across frames while the same features touch
};

struct ContactManifold {
    Vec3 normal;  // unit, points from shape A toward shape B
    ManifoldPoint points[8];
    int count = 0;
};

// Builds the manifold for shapes A and B. Returns false when they are separated by more than `margin`.
bool collide(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb, float margin,
             ContactManifold& out);

AABB computeAABB(const Shape& s, const Transform& t);

}  // namespace rbe
