#include "PairOrderer.h"
#include "BodySystem.h"
#include "debug/debug.h"
#include "Body.h"

#include "render/Graphics.h"

#include <vector>
#include <algorithm>

extern Graphics* g_graphics;

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
        return sum_mass1 > sum_mass2;
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
        return sum_mass1 > sum_mass2;
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
    // do_not_order();

    order_by_mass_and_tree();
    // tree.print();

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
#if 0
        DebugBegin();
        dout << green << "After sorting transitions.." << newl;
        for (auto x : transition_list.second) {
            Body b1 = body_system.get_body( broadphase_alg.get_box_user_data(tree.nodes[x].first) );
            Body b2 = body_system.get_body( broadphase_alg.get_box_user_data(tree.nodes[x].second) );
            dout << b1.shape().get_mass() + b2.shape().get_mass() << newl;
        }
#endif
    }
    traverse();
}


// Depth-first //
void PairOrderer::traverse()
{
    DebugBeginC(false);
    traversed_pairs.clear();

    bool* visited = new bool[tree.nodes.size()] ();     // nodes visited

    int current_node;

    for ( int indices_index = 0; indices_index < sorted_indices.size(); ++ indices_index )
    {
        dout << "Back to main loop.." << newl;
        current_node = sorted_indices[indices_index];
        if ( !visited[current_node] )
            traverse( current_node, visited );
    }
    delete[] visited;
}

void PairOrderer::traverse( int node, bool* visited )
{
    DebugBegin();
    visited[node] = true;

    dout << "Visit node " << tree.nodes[node].first << ", " << tree.nodes[node].second << newl;       
    traversed_pairs.push_back( tree.nodes[node] );
    traverse_neighbors( node, visited );
}

void PairOrderer::traverse_neighbors( int node, bool* visited )
{
    // Traverse the neighbors of the box with the least mass first
    if (body_system.get_body(tree.nodes[node].first).shape().get_mass()
            < body_system.get_body(tree.nodes[node].second).shape().get_mass()) {
        traverse_neighbors_of_box(node, tree.nodes[node].first, visited);
        traverse_neighbors_of_box(node, tree.nodes[node].second, visited);
    } else {
        traverse_neighbors_of_box(node, tree.nodes[node].second, visited);
        traverse_neighbors_of_box(node, tree.nodes[node].first, visited);
    }
}
void PairOrderer::traverse_neighbors_of_box( int node, int box_id, bool* visited)
{
    if (tree.used_in_nodes.count(   box_id  )) {                  // if there are adjacent pairs
        for (int i : tree.used_in_nodes[    box_id  ]) {          // loop through them & add to open set
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
