#include "PairOrderer.h"
#include "BodySystem.h"
#include "debug/debug.h"
#include "Body.h"

#include "render/Graphics.h"

#include <vector>
#include <algorithm>

extern Renderer* g_graphics;

// struct for ordering by sum of masses of a pair of bodies
struct lt_indices_by_mass
{
 public:
    lt_indices_by_mass(PairTree& tree, BodySystem& body_system)
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
struct lt_pairs_by_mass
{
    lt_pairs_by_mass(BodySystem& body_system) : body_system(body_system) {}
    bool operator()(Pair pair1, Pair pair2) {
        const float sum_mass1 = body_system.get_body(pair1.first).shape().get_mass()
                              + body_system.get_body(pair1.second).shape().get_mass();
        const float sum_mass2 = body_system.get_body(pair2.first).shape().get_mass()
                              + body_system.get_body(pair2.second).shape().get_mass();
        return sum_mass1 < sum_mass2;
    }
 private:
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
    do_not_order();

    // order_by_mass_and_tree();

    // not effective at all:
    // simple_order_by_mass();
}

void PairOrderer::do_not_order()
{
    traversed_pairs.clear();
    // TODO (performance) resize capacity
    for (Pair pair : broadphase_alg.pairs) {
        traversed_pairs.push_back(pair);
    }
}
void PairOrderer::simple_order_by_mass()
{
    traversed_pairs.clear();
    // TODO (performance) resize capacity
    for (Pair pair : broadphase_alg.pairs) {
        traversed_pairs.push_back(pair);
    }
    lt_pairs_by_mass comparator(body_system);
    std::sort(traversed_pairs.begin(), traversed_pairs.end(), comparator);
}

std::vector<int> create_range(int size)
{
    std::vector<int> result;
    for (int i = 0; i < size; i ++) {
        result.push_back(i);
    }
    return result;
}

void PairOrderer::order_by_mass_and_tree()
{
    tree = {};
    for (Pair pair : broadphase_alg.pairs) {
        // Convert from box index to body index //
        pair.first = broadphase_alg.get_box_user_data(pair.first);
        pair.second = broadphase_alg.get_box_user_data(pair.second);
        tree.add_pair(pair);
    }


    // Sort the tree's node indices by mass //
    lt_indices_by_mass comparator(tree, body_system);
    sorted_indices = create_range(tree.num_nodes());
    std::sort(sorted_indices.begin(), sorted_indices.end(), comparator);

    // Also sort the transitions/edges by mass //
    for (auto transition_list : tree.used_in_nodes)
    {
        std::sort(transition_list.second.begin(), transition_list.second.end(), comparator);
    }
    traverse();
}


// Depth-first //
void PairOrderer::traverse()
{
    DebugBegin();
    traversed_pairs.clear();

    bool* visited = new bool[tree.nodes.size()] ();     // nodes visited

    int current_node;

    for ( int indices_index = 0; indices_index < sorted_indices.size(); ++ indices_index )
    {
        current_node = sorted_indices[indices_index];
        if ( !visited[current_node] )
            traverse( current_node, visited );
    }
    delete[] visited;
}

void PairOrderer::traverse( int node, bool* visited )
{
    visited[node] = true;

    traversed_pairs.push_back( tree.nodes[node] );
    traverse_neighbors( node, visited );
}

void PairOrderer::traverse_neighbors( int node, bool* visited )
{
    if (tree.used_in_nodes.count(   tree.nodes[node].first  )) {                // if there are adjacent pairs
        for (int i : tree.used_in_nodes[    tree.nodes[node].first  ]) {        // loop through them & add to open set
            if (!visited[i] && i != node)
                traverse(i, visited);
        }
    }
    if (tree.used_in_nodes.count(   tree.nodes[node].second  )) {
        for (int i : tree.used_in_nodes[    tree.nodes[node].second  ]) {
            if (!visited[i] && i != node)
                traverse(i, visited);
        }
    }
}


void PairOrderer::visualize()
{
    for (Pair pair : traversed_pairs)
    {
    }
}
