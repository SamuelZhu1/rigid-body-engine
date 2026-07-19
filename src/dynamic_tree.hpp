// Dynamic AABB tree broadphase (same idea as Box2D's b2DynamicTree).
// Leaves hold fat AABBs so a body that jiggles in place doesn't get re-inserted every step.
#pragma once

#include <cstdint>
#include <vector>

#include "rbe/math.hpp"

namespace rbe {

class DynamicTree {
public:
    static constexpr int32_t kNull = -1;

    int32_t createProxy(const AABB& tight, float margin, uint32_t userData);
    void destroyProxy(int32_t proxy);
    // Returns true if the proxy was re-inserted (tight box left the fat box).
    bool moveProxy(int32_t proxy, const AABB& tight, float margin);

    const AABB& fatAABB(int32_t proxy) const { return nodes_[proxy].box; }
    uint32_t userData(int32_t proxy) const { return nodes_[proxy].userData; }

    // Calls f(userData) for every leaf whose fat box overlaps `box`.
    template <typename F>
    void query(const AABB& box, F&& f) const {
        if (root_ == kNull) return;
        stack_.clear();
        stack_.push_back(root_);
        while (!stack_.empty()) {
            int32_t id = stack_.back();
            stack_.pop_back();
            const Node& n = nodes_[id];
            if (!n.box.overlaps(box)) continue;
            if (n.isLeaf()) {
                f(n.userData);
            } else {
                stack_.push_back(n.child1);
                stack_.push_back(n.child2);
            }
        }
    }

    int height() const { return root_ == kNull ? 0 : nodes_[root_].height; }
    int leafCount() const { return leafCount_; }
    // Recomputes heights and boxes from scratch and checks parent links; used by tests.
    bool validate() const;
    void collectNodeBoxes(std::vector<AABB>& out) const;

private:
    struct Node {
        AABB box;
        int32_t parent = kNull;  // doubles as the free-list "next" pointer
        int32_t child1 = kNull, child2 = kNull;
        int32_t height = -1;  // leaf = 0, free = -1
        uint32_t userData = 0;
        bool isLeaf() const { return child1 == kNull; }
    };

    int32_t allocateNode();
    void freeNode(int32_t id);
    void insertLeaf(int32_t leaf);
    void removeLeaf(int32_t leaf);
    int32_t balance(int32_t a);
    int validateNode(int32_t id, bool& ok) const;

    std::vector<Node> nodes_;
    int32_t root_ = kNull;
    int32_t freeList_ = kNull;
    int leafCount_ = 0;
    mutable std::vector<int32_t> stack_;
};

}  // namespace rbe
