#include "PairTree.h"
#include "PairManager.h" // for Pair

void PairTree::add_pair(Pair pair)
{
    // TODO: could update neighbors, if we intend to use it. Maybe better to just use used_in_nodes

    int index = nodes.size();
    nodes.push_back(pair);
    used_in_nodes[pair.first].push_back(index);
    used_in_nodes[pair.second].push_back(index);
}
