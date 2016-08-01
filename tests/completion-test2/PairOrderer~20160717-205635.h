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
    std::vector<Pair>::iterator begin();
    std::vector<Pair>::iterator end();
 private:
    void order();
    void traverse();
    void add_neighbors_to_open_set(int node, bool* visited, std::deque<int>& open_set);
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
