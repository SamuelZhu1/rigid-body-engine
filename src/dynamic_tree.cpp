#include "dynamic_tree.hpp"

#include <algorithm>
#include <cassert>

namespace rbe {

int32_t DynamicTree::allocateNode() {
    if (freeList_ == kNull) {
        nodes_.emplace_back();
        return static_cast<int32_t>(nodes_.size() - 1);
    }
    int32_t id = freeList_;
    freeList_ = nodes_[id].parent;
    nodes_[id] = Node{};
    return id;
}

void DynamicTree::freeNode(int32_t id) {
    nodes_[id].parent = freeList_;
    nodes_[id].height = -1;
    freeList_ = id;
}

int32_t DynamicTree::createProxy(const AABB& tight, float margin, uint32_t userData) {
    int32_t id = allocateNode();
    nodes_[id].box = tight.fattened(margin);
    nodes_[id].userData = userData;
    nodes_[id].height = 0;
    insertLeaf(id);
    ++leafCount_;
    return id;
}

void DynamicTree::destroyProxy(int32_t proxy) {
    removeLeaf(proxy);
    freeNode(proxy);
    --leafCount_;
}

bool DynamicTree::moveProxy(int32_t proxy, const AABB& tight, float margin) {
    if (nodes_[proxy].box.contains(tight)) return false;
    removeLeaf(proxy);
    nodes_[proxy].box = tight.fattened(margin);
    insertLeaf(proxy);
    return true;
}

void DynamicTree::insertLeaf(int32_t leaf) {
    if (root_ == kNull) {
        root_ = leaf;
        nodes_[leaf].parent = kNull;
        return;
    }

    // Walk down picking the child that adds the least surface area (branch and bound lite).
    const AABB leafBox = nodes_[leaf].box;
    int32_t index = root_;
    while (!nodes_[index].isLeaf()) {
        const Node& n = nodes_[index];
        float area = n.box.surfaceArea();
        float combined = merge(n.box, leafBox).surfaceArea();
        float costHere = 2.0f * combined;             // make a new parent here
        float inherit = 2.0f * (combined - area);     // cost pushed down to children

        auto childCost = [&](int32_t c) {
            const Node& cn = nodes_[c];
            float m = merge(leafBox, cn.box).surfaceArea();
            return cn.isLeaf() ? m + inherit : (m - cn.box.surfaceArea()) + inherit;
        };
        float c1 = childCost(n.child1), c2 = childCost(n.child2);
        if (costHere < c1 && costHere < c2) break;
        index = c1 < c2 ? n.child1 : n.child2;
    }

    int32_t sibling = index;
    int32_t oldParent = nodes_[sibling].parent;
    int32_t newParent = allocateNode();
    nodes_[newParent].parent = oldParent;
    nodes_[newParent].box = merge(leafBox, nodes_[sibling].box);
    nodes_[newParent].height = nodes_[sibling].height + 1;
    nodes_[newParent].child1 = sibling;
    nodes_[newParent].child2 = leaf;
    nodes_[sibling].parent = newParent;
    nodes_[leaf].parent = newParent;

    if (oldParent != kNull) {
        if (nodes_[oldParent].child1 == sibling)
            nodes_[oldParent].child1 = newParent;
        else
            nodes_[oldParent].child2 = newParent;
    } else {
        root_ = newParent;
    }

    // Refit and rebalance up to the root.
    index = nodes_[leaf].parent;
    while (index != kNull) {
        index = balance(index);
        Node& n = nodes_[index];
        n.height = 1 + std::max(nodes_[n.child1].height, nodes_[n.child2].height);
        n.box = merge(nodes_[n.child1].box, nodes_[n.child2].box);
        index = n.parent;
    }
}

void DynamicTree::removeLeaf(int32_t leaf) {
    if (leaf == root_) {
        root_ = kNull;
        return;
    }
    int32_t parent = nodes_[leaf].parent;
    int32_t grand = nodes_[parent].parent;
    int32_t sibling = nodes_[parent].child1 == leaf ? nodes_[parent].child2 : nodes_[parent].child1;

    if (grand != kNull) {
        if (nodes_[grand].child1 == parent)
            nodes_[grand].child1 = sibling;
        else
            nodes_[grand].child2 = sibling;
        nodes_[sibling].parent = grand;
        freeNode(parent);

        int32_t index = grand;
        while (index != kNull) {
            index = balance(index);
            Node& n = nodes_[index];
            n.box = merge(nodes_[n.child1].box, nodes_[n.child2].box);
            n.height = 1 + std::max(nodes_[n.child1].height, nodes_[n.child2].height);
            index = n.parent;
        }
    } else {
        root_ = sibling;
        nodes_[sibling].parent = kNull;
        freeNode(parent);
    }
}

// AVL-style rotation: if one child is 2+ levels taller, promote its taller grandchild.
int32_t DynamicTree::balance(int32_t iA) {
    Node* A = &nodes_[iA];
    if (A->isLeaf() || A->height < 2) return iA;

    int32_t iB = A->child1, iC = A->child2;
    int32_t diff = nodes_[iC].height - nodes_[iB].height;

    auto rotate = [&](int32_t iUp, int32_t iOther) {
        // iUp is the taller child of A; its taller child stays, the other moves under A.
        Node& U = nodes_[iUp];
        int32_t iF = U.child1, iG = U.child2;
        U.child1 = iA;
        U.parent = nodes_[iA].parent;
        nodes_[iA].parent = iUp;

        if (U.parent != kNull) {
            if (nodes_[U.parent].child1 == iA)
                nodes_[U.parent].child1 = iUp;
            else
                nodes_[U.parent].child2 = iUp;
        } else {
            root_ = iUp;
        }

        Node& Aref = nodes_[iA];
        int32_t keep = nodes_[iF].height > nodes_[iG].height ? iF : iG;
        int32_t move = keep == iF ? iG : iF;
        U.child2 = keep;
        if (Aref.child1 == iUp)
            Aref.child1 = move;
        else
            Aref.child2 = move;
        nodes_[move].parent = iA;
        Aref.box = merge(nodes_[iOther].box, nodes_[move].box);
        Aref.height = 1 + std::max(nodes_[iOther].height, nodes_[move].height);
        U.box = merge(Aref.box, nodes_[keep].box);
        U.height = 1 + std::max(Aref.height, nodes_[keep].height);
        return iUp;
    };

    if (diff > 1) return rotate(iC, iB);
    if (diff < -1) return rotate(iB, iC);
    return iA;
}

int DynamicTree::validateNode(int32_t id, bool& ok) const {
    const Node& n = nodes_[id];
    if (n.isLeaf()) {
        if (n.height != 0) ok = false;
        return 0;
    }
    if (nodes_[n.child1].parent != id || nodes_[n.child2].parent != id) ok = false;
    int h1 = validateNode(n.child1, ok), h2 = validateNode(n.child2, ok);
    int h = 1 + std::max(h1, h2);
    if (h != n.height) ok = false;
    if (!n.box.contains(nodes_[n.child1].box) || !n.box.contains(nodes_[n.child2].box)) ok = false;
    return h;
}

bool DynamicTree::validate() const {
    if (root_ == kNull) return leafCount_ == 0;
    bool ok = nodes_[root_].parent == kNull;
    validateNode(root_, ok);
    return ok;
}

void DynamicTree::collectNodeBoxes(std::vector<AABB>& out) const {
    if (root_ == kNull) return;
    std::vector<int32_t> st{root_};
    while (!st.empty()) {
        int32_t id = st.back();
        st.pop_back();
        out.push_back(nodes_[id].box);
        if (!nodes_[id].isLeaf()) {
            st.push_back(nodes_[id].child1);
            st.push_back(nodes_[id].child2);
        }
    }
}

}  // namespace rbe
