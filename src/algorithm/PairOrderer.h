#pragma once
#include "SAP.h"
#include "PairTree.h"
#include <vector>
#include <deque>

class BodySystem;

/**
For efficiency, orders the pairs that the broadphase algorithm yields.
Orders by mass.
**/

class PairOrderer
{
 public:
    PairOrderer(SAP<int>& broadphase_alg, BodySystem& body_system);
    void update();
    void visualize();
    std::vector<Pair>::iterator begin();
    std::vector<Pair>::iterator end();
 private:
    // Order algorithms //
    void do_not_order();
    void simple_order_by_mass();
    void order_by_mass_and_tree();

    // Helpers //
    void traverse();
    void traverse(int node, bool* visited);
    void traverse_neighbors(int node, bool* visited);
 private:
    // Persistent:
    SAP<int>& broadphase_alg;
    BodySystem& body_system;
    // After ordering:
    PairTree tree;
    std::vector<int> sorted_indices;
    // After traversing:
    std::vector<Pair> traversed_pairs;

};
