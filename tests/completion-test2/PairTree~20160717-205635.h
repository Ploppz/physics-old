#pragma once
#include "PairManager.h"
#include <vector>
#include <unordered_map>

class PairOrderer;

// Note: `id` in the comments refers to an element/member of a pair

class PairTree
{
 public:
    void add_pair(Pair pair);

    Pair get_node(int index) { return nodes[index]; }
    int num_nodes()         { return nodes.size(); }

    friend PairOrderer;
 private:

    // translates from node to index
    /* std::unordered_map<Pair, int> index_map; */

    // translates from index to node
    std::vector<Pair> nodes;

    // [id] -> [list of indices of nodes containing this id]
    std::unordered_map<int, std::vector<int>> used_in_nodes;
};
