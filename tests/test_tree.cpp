#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "dynamic_tree.hpp"
#include "rbe/world.hpp"
#include "scenes.hpp"

using namespace rbe;

namespace {
struct Lcg {
    uint32_t s;
    float next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<float>(s >> 8) / 16777216.0f;
    }
};
AABB randomBox(Lcg& r, float world) {
    Vec3 c{(r.next() - 0.5f) * world, (r.next() - 0.5f) * world, (r.next() - 0.5f) * world};
    Vec3 h{0.2f + r.next(), 0.2f + r.next(), 0.2f + r.next()};
    return {c - h, c + h};
}
}  // namespace

TEST(DynamicTree, QueryMatchesBruteForceAfterRandomMoves) {
    Lcg rng{7};
    DynamicTree tree;
    const int n = 1000;
    std::vector<AABB> tight(n);
    std::vector<int32_t> proxies(n);
    for (int i = 0; i < n; ++i) {
        tight[i] = randomBox(rng, 60);
        proxies[i] = tree.createProxy(tight[i], 0.1f, static_cast<uint32_t>(i));
    }
    for (int round = 0; round < 100; ++round) {
        for (int k = 0; k < 50; ++k) {
            int i = static_cast<int>(rng.next() * n);
            Vec3 d{(rng.next() - 0.5f) * 2, (rng.next() - 0.5f) * 2, (rng.next() - 0.5f) * 2};
            tight[i] = {tight[i].lo + d, tight[i].hi + d};
            tree.moveProxy(proxies[i], tight[i], 0.1f);
        }
        // Every fat-box overlap the tree reports must equal the brute-force set.
        std::set<std::pair<int, int>> fromTree, brute;
        for (int i = 0; i < n; ++i) {
            tree.query(tree.fatAABB(proxies[i]), [&](uint32_t j) {
                if (static_cast<int>(j) > i) fromTree.emplace(i, static_cast<int>(j));
            });
            for (int j = i + 1; j < n; ++j)
                if (tree.fatAABB(proxies[i]).overlaps(tree.fatAABB(proxies[j]))) brute.emplace(i, j);
        }
        ASSERT_EQ(fromTree, brute) << "round " << round;
    }
    EXPECT_TRUE(tree.validate());
}

TEST(DynamicTree, StaysBalanced) {
    Lcg rng{3};
    DynamicTree tree;
    const int n = 4096;
    // Inserting boxes along a line is the worst case for an unbalanced tree.
    for (int i = 0; i < n; ++i) {
        AABB b{{i * 1.0f, 0, 0}, {i * 1.0f + 0.5f, 0.5f, 0.5f}};
        tree.createProxy(b, 0.0f, static_cast<uint32_t>(i));
    }
    EXPECT_TRUE(tree.validate());
    EXPECT_LE(tree.height(), 2 * static_cast<int>(std::log2(n)) + 1);
}

TEST(DynamicTree, DestroyKeepsTreeValid) {
    Lcg rng{11};
    DynamicTree tree;
    std::vector<int32_t> ids;
    for (int i = 0; i < 300; ++i) ids.push_back(tree.createProxy(randomBox(rng, 30), 0.1f, static_cast<uint32_t>(i)));
    for (int i = 0; i < 300; i += 2) tree.destroyProxy(ids[i]);
    EXPECT_EQ(tree.leafCount(), 150);
    EXPECT_TRUE(tree.validate());
}

TEST(Broadphase, TreeAndBruteForceSimulateIdentically) {
    // Same pairs in the same order means bit-identical simulations.
    auto run = [](BroadphaseMode mode) {
        WorldSettings s;
        s.broadphase = mode;
        World w(s);
        scenes::build(w, "pile", 300);
        for (int i = 0; i < 200; ++i) w.step(1.0f / 60.0f);
        return w.stateHash();
    };
    EXPECT_EQ(run(BroadphaseMode::Tree), run(BroadphaseMode::BruteForce));
}
