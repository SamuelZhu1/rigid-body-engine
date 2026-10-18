#include "scenes.hpp"

#include <cmath>
#include <cstdint>

namespace rbe::scenes {

namespace {
// Tiny deterministic generator so scenes are identical on every platform (std:: distributions are not).
struct Lcg {
    uint32_t s;
    float next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<float>(s >> 8) / 16777216.0f;  // [0, 1)
    }
};

BodyId staticBox(World& w, const Vec3& p, const Vec3& half) {
    BodyDef d;
    d.shape.type = ShapeType::Box;
    d.shape.halfExtents = half;
    d.position = p;
    d.mass = 0;
    return w.createBody(d);
}
}  // namespace

BodyId addGround(World& w) { return staticBox(w, {0, -0.5f, 0}, {100, 0.5f, 100}); }

std::vector<BodyId> stack(World& w, int n) {
    std::vector<BodyId> ids;
    for (int i = 0; i < n; ++i) {
        BodyDef d;
        d.position = {0, 0.5f + 1.0f * i, 0};
        ids.push_back(w.createBody(d));
    }
    return ids;
}

std::vector<BodyId> pyramid(World& w, int rows) {
    std::vector<BodyId> ids;
    for (int r = 0; r < rows; ++r) {
        int count = rows - r;
        for (int i = 0; i < count; ++i) {
            BodyDef d;
            d.position = {(i - (count - 1) * 0.5f) * 1.05f, 0.5f + r * 1.0f, 0};
            ids.push_back(w.createBody(d));
        }
    }
    return ids;
}

std::vector<BodyId> pile(World& w, int n, unsigned seed) {
    // Pit sized so the boxes settle about 6 layers deep.
    int side = static_cast<int>(std::ceil(std::sqrt(n / 6.0f)));
    float half = side * 0.6f + 0.5f;
    staticBox(w, {half + 0.5f, 5, 0}, {0.5f, 5, half + 1});
    staticBox(w, {-half - 0.5f, 5, 0}, {0.5f, 5, half + 1});
    staticBox(w, {0, 5, half + 0.5f}, {half + 1, 5, 0.5f});
    staticBox(w, {0, 5, -half - 0.5f}, {half + 1, 5, 0.5f});

    Lcg rng{seed};
    std::vector<BodyId> ids;
    for (int i = 0; i < n; ++i) {
        int layer = i / (side * side), k = i % (side * side);
        BodyDef d;
        d.position = {(k % side - (side - 1) * 0.5f) * 1.2f, 1.0f + layer * 1.3f, (k / side - (side - 1) * 0.5f) * 1.2f};
        d.orientation = Quat::fromAxisAngle({rng.next() - 0.5f, rng.next() - 0.5f, rng.next() - 0.5f}, rng.next() * 3.0f);
        ids.push_back(w.createBody(d));
    }
    return ids;
}

BodyId spinningBox(World& w, float spin) {
    w.settings().gravity = {0, 0, 0};
    BodyDef d;
    d.shape.halfExtents = {1.0f, 0.5f, 0.2f};  // three different moments of inertia
    d.position = {0, 3, 0};
    // Mostly about the middle (y) axis, plus a small nudge that lets the instability grow.
    d.angularVelocity = {0.01f * spin, spin, 0.01f * spin};
    return w.createBody(d);
}

std::vector<BodyId> dominoes(World& w, int n) {
    std::vector<BodyId> ids;
    for (int i = 0; i < n; ++i) {
        BodyDef d;
        d.shape.halfExtents = {0.1f, 1.0f, 0.5f};
        d.position = {i * 1.2f, 1.0f, 0};
        d.friction = 0.4f;
        ids.push_back(w.createBody(d));
    }
    // Tip the first one.
    w.setVelocity(ids[0], {0, 0, 0}, {0, 0, -2.0f});
    return ids;
}

std::vector<BodyId> sphereRain(World& w, int n, unsigned seed) {
    Lcg rng{seed};
    std::vector<BodyId> ids;
    for (int i = 0; i < n; ++i) {
        BodyDef d;
        d.shape.type = ShapeType::Sphere;
        d.shape.radius = 0.3f + 0.3f * rng.next();
        d.position = {(rng.next() - 0.5f) * 10, 2 + i * 0.5f, (rng.next() - 0.5f) * 10};
        ids.push_back(w.createBody(d));
    }
    return ids;
}

bool build(World& w, const std::string& name, int count) {
    if (name == "spin") {
        spinningBox(w);
        return true;
    }
    addGround(w);
    if (name == "stack") stack(w, count > 0 ? count : 20);
    else if (name == "pyramid") pyramid(w, count > 0 ? count : 20);
    else if (name == "pile") pile(w, count > 0 ? count : 2000);
    else if (name == "dominoes") dominoes(w, count > 0 ? count : 30);
    else if (name == "spheres") sphereRain(w, count > 0 ? count : 500);
    else return false;
    return true;
}

}  // namespace rbe::scenes
