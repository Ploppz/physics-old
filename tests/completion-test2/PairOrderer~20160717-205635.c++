#include "PairOrderer.h"
// #include "BodySystem.h"
// #include "debug/debug.h"
// #include "Body.h"
#include <vector>
#include <algorithm>
// struct for ordering by sum of masses of a pair of bodies
struct lt_by_mass
{
 public:
    lt_by_mass(PairTree& tree, BodySystem& body_system)
        : tree(tree), body_system(body_system) {}
    bool operator()(int index1, int index2) {
        const float sum_mass1 = body_system.get_body(tree.get_node(index1).first).shape().get_mass()
                              + body_system.get_body(tree.get_node(index1).second).shape().get_mass();
        const float sum_mass2 = body_system.get_body(tree.get_node(index2).first).shape().get_mass()
                              + body_system.get_body(tree.get_node(index2).second).shape().get_mass();
        return sum_mass1 < sum_mass2;
    }
 private:
    PairTree& tree;
    BodySystem& body_system;
};


PairOrderer::PairOrderer(SAP<int>& broadphase_alg, BodySystem& body_system)
    : broadphase_alg(broadphase_alg), body_system(body_system)
{
}


std::vector<Pair>::iterator PairOrderer::begin()
{
    return traversed_pairs.begin();
}
std::vector<Pair>::iterator PairOrderer::end()
{
    return traversed_pairs.end();
}


void PairOrderer::update()
{
    order();
    traverse();
}

std::vector<int> create_range(int size)
{
    std::vector<int> result;
    for (int i = 0; i < size; i ++) {
        result.push_back(i);
    }
    return result;
}

void PairOrderer::order()
{
    tree = {};
    for (Pair pair : broadphase_alg.pairs) {
        // Convert from box index to body index //
        pair.first = broadphase_alg.get_box_user_data(pair.first);
        pair.second = broadphase_alg.get_box_user_data(pair.second);
        tree.add_pair(pair);
    }

    lt_by_mass comparator(tree, body_system);

    sorted_indices = create_range(tree.num_nodes());
    std::sort(sorted_indices.begin(), sorted_indices.end(), comparator);

}


// Breadth-first //
void PairOrderer::traverse()
{
    DebugBegin();
    traversed_pairs.clear();
    bool* visited = new bool[tree.nodes.size()] ();     // nodes visited
    std::deque<int> open_set;      // indices into the tree nodes
    int current_node;

    for (int indices_index = 0; indices_index < sorted_indices.size(); ++ indices_index)
    {
        current_node = sorted_indices[indices_index];
        if (visited[current_node]) continue;
        visited[current_node] = true;

        dout << "current node: " << current_node << newl;
        add_neighbors_to_open_set(current_node, visited, open_set);
        traversed_pairs.push_back(tree.nodes[ current_node ]);

        while (open_set.size() > 0)
        {
            current_node = open_set.front(); open_set.pop_front();
            if (visited[current_node]) continue;
            visited[current_node] = true;

            add_neighbors_to_open_set(current_node, visited, open_set);
            traversed_pairs.push_back(tree.nodes[ current_node ]);
        }
    }
    delete[] visited;
}

void PairOrderer::add_neighbors_to_open_set(int node, bool* visited, std::deque<int>& open_set)
{
    if (tree.used_in_nodes.count(   tree.nodes[node].first  )) {                // if there are adjacent pairs
        for (int i : tree.used_in_nodes[    tree.nodes[node].first  ]) {        // loop through them & add to open set
            if (!visited[i] && i != node)
                open_set.push_back(i);
        }
    }
    if (tree.used_in_nodes.count(   tree.nodes[node].second  )) {
        for (int i : tree.used_in_nodes[    tree.nodes[node].second  ]) {
            if (!visited[i] && i != node)
                open_set.push_back(i);
        }
    }
}
