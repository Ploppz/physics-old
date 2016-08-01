#include "PairTree.h"
#include "PairManager.h" // for Pair
#include "debug/debug.h"

void PairTree::add_pair(Pair pair)
{
    // TODO: could update neighbors, if we intend to use it. Maybe better to just use used_in_nodes

    int index = nodes.size();
    nodes.push_back(pair);
    used_in_nodes[pair.first].push_back(index);
    used_in_nodes[pair.second].push_back(index);
}

void PairTree::print()
{
    DebugBegin();
    dout << "Tree:" << newl;
    dout << "  Pairs: " << newl;
    int counter = 0;
    for (Pair p : nodes)
    {
        dout << "   #" << counter << ": " << p.first << ", " << p.second << newl;
        ++ counter;
    }
    dout << "  Transitions: " << newl;
    for (auto p : used_in_nodes) {
        dout << "   box " << p.first << " -> pairs: ";
        for (auto it : p.second)
        {
            dout << it << ",";
        }
        dout << newl;
    }
}
