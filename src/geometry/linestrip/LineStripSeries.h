#pragma once
#include "LineStrip.h"
#include "linestrip_iterator.h"
#include "../geometry.h"
#include <list>
#include <glm/glm.hpp>
#include <iostream>


class Polygon;

template<Direction x_bias>
class LineStripSeries
{
 public:
    template <bool reverse>
    class Vertex;
    //
    typedef std::list<LineStrip>::iterator iterator;
    LineStripSeries(LineStrip to_copy);

    //
    void make_y_monotone();
    void append_lines_to_vector(std::vector<float> &list);
    void set_align_matrix(glm::mat2 matrix) { align_transform = matrix; };
    Polygon* get_parent() { return parent; }

    template <bool reverse>
    void print_with_current(Vertex<reverse> curr);


 private:
    template <bool reverse>
    void make_y_monotone_direction();

    ////
    glm::mat2 align_transform;
    Polygon* parent;
 public:
    std::list<LineStrip> linestrips;
};


#define LSS LineStripSeries<x_bias>

template <Direction x_bias>
template <bool reverse>
class LSS::Vertex
{
 public:
    /* Vertex() {} */
    Vertex(LineStripSeries* parent);

    EdgePoint to_edge_point(float alpha)
    {
        return vertex.to_edge_point(alpha);
    }

    /* Information */
    // TODO For now, all these are in transformed and aligned space (uses parent->align_transform)
    
    glm::vec2 entry_point();
    glm::vec2 exit_point();

    // Picks one of entry or exit points based on their x values and the x_bias
    glm::vec2 best_point();

    bool at_end();
    bool at_start();

    /* Operations */
    void set_to_start();
    LSS::Vertex<reverse>& operator ++();
    LSS::Vertex<reverse>& operator --();

    void split_linestrip_and_move_to_new();

    // Set the start of the current linestrip to a place between this vertex and the next given by alpha [0, 1]
    void set_linestrip_start(float alpha);

    LineStrip::Vertex<reverse> get_vertex() { return vertex; }
    linestrip_iterator<reverse> get_linestrip_it() { return linestrip_it; }
    LineStripSeries* get_parent() { return parent; };
  /* private: */
    LineStripSeries* parent;
    linestrip_iterator<reverse> linestrip_it;
    LineStrip::Vertex<reverse> vertex;
};



/***** LineStripSeries IMPLEMENTATION *****/

template <Direction x_bias>
std::ostream& operator<< (std::ostream& out, LineStripSeries<x_bias>& lss)
{
    int c = 0;
    for (auto it = lss.linestrips.begin(); it != lss.linestrips.end(); ++ it)
    {
        c++;
        out << "\t" << *it << std::endl;
    }
    return out;
}

template<Direction x_bias>
LSS::LineStripSeries(LineStrip to_copy)
{
    this->parent = to_copy.get_parent();
    linestrips.push_back(to_copy);
}
template <Direction x_bias>
void LineStripSeries<x_bias>::make_y_monotone()
{
    std::cout << "\n\n" << std::endl;
    std::cout << "### Input: " << *this << std::endl;
    std::cout << "### Monotonizing Forward ###" << std::endl;
    make_y_monotone_direction<false>();
    std::cout << "### Monotonizing Backward ###" << std::endl;
    make_y_monotone_direction<true>();
    
}

/***
In order to make it clean, I don't think about performance.
We have 4 vertices (2 which may be the same) and 2 edges to keep track of.
It's possible to 'cache' these as looping variables.
Takes more space, but avoids a lot of transformations.
***/
template<Direction x_bias>
template<bool reverse>
void LineStripSeries<x_bias>::make_y_monotone_direction()
{
    /* Initiate */
    Vertex<reverse> curr(this);
    curr.set_to_start();

    Vertex<reverse> next = curr; ++ next;
    glm::vec2 a = next.entry_point() - curr.exit_point();

    curr = next; ++ next;
    glm::vec2 b = next.entry_point() - curr.exit_point();
    /**/

    bool covered = false;
    float uncover_y = 0;

    int count = 0;
    while (true) {
        if (next.at_end())
            break;
        ++ count;
        if (count > 200) {
            exit(0); 
        }
        /* print_with_current(curr); // debug  */
        if ( ! covered) {
            bool b_prone_to_be_hidden = leftof(b, a) ^ (a.y < 0) ^ (x_bias == LEFT);    // leftof(a, b)   if   a.y < 0 XOR x == LEFT
            bool b_and_a_opposite_directions = (b.y < 0) ^ (a.y < 0);                   // b.y >= 0       if   a.y < 0
            
            if (b_prone_to_be_hidden && b_and_a_opposite_directions) { // Edge b is hidden.
                uncover_y = curr.exit_point().y;
                curr.split_linestrip_and_move_to_new();
                /* Update next... */
                next = curr; ++ next;
                covered = true;
                /* marks the start of deletion */
            }
        } else {
            float t;
            intersect_horizontal(curr.exit_point(), b, uncover_y, t);
            if (t >= 0 && t < 1) {
                // Edge b is partially or wholly uncovered. This is where the new LS can start.
                curr.set_linestrip_start(t);
                covered = false;
                /* marks the end of deletion */
            }
        }
        // Advance
        if (next.at_end())
            break;
        curr = next; ++ next;
        a = b;
        b = next.entry_point() - curr.exit_point();
    }
}

template <Direction x_bias>
void LineStripSeries<x_bias>::append_lines_to_vector(std::vector<float> &list)
{
    /* std::cout << "RENDER ** LineStripSeries" << std::endl; */
    Vertex<false> curr(this);
    curr.set_to_start();
    Vertex<false> next = curr;
    ++ next;
    while (true) {
        /* std::cout << "\t at " << *(curr.linestrip_it.iterator) << std::endl;
        std::cout << "\t\t" << curr.vertex.get_index() << std::endl;
        std::cout << ".. " << curr.exit_point() << ", " << next.entry_point() << std::endl; */
        glm::vec2 vec_i = curr.exit_point();
        glm::vec2 vec_j = next.entry_point();

        list.push_back(vec_i.x);
        list.push_back(vec_i.y);

        list.push_back(vec_j.x);
        list.push_back(vec_j.y);

        if (next.at_end()) break;
        curr = next;
        ++ next;
    }
    /* std::cout << "END RENDER ** LineStripSeries" << std::endl; */
}



template <Direction x_bias>
template <bool reverse>
void LineStripSeries<x_bias>::print_with_current(Vertex<reverse> curr)
{
    int c = 0;
    std::cout << "LSS {" << std::endl;
    for (auto it = linestrips.begin(); it != linestrips.end(); it ++)
    {
        c++;
        if (*it == *curr.linestrip_it.base()) {
            std::cout << "\t**" << *it << std::endl;
        } else {
            std::cout << "\t" << *it << std::endl;
        }
    }
    std::cout << "}" << std::endl;
}

/////////////////////////////////
/** Vertex of LineStripSeries **/
/////////////////////////////////

/**
 The Vertex is never at the end of a LineStrip, except, on the
very last LineStrip in the list.
If it's at the start of a LineStrip (except the first one),
then entry_point will return the end of the previous LineStrip,
and exit_point will return the start of the current one.

Another note:
 All references to a 'start' or an 'end' takes into account whether
reverse == true.
**/

#define templates template <Direction x_bias> template <bool reverse>

templates
LSS::Vertex<reverse>::Vertex(LineStripSeries* parent)
    : parent(parent)
{
    linestrip_it.set_to_start_of(parent->linestrips);
    // vertex.parent = &(*linestrip_it.iterator);
    vertex.parent = &(* parent->linestrips.begin());
    vertex.set_to_start();
}
templates
void LSS::Vertex<reverse>::split_linestrip_and_move_to_new()
{
    // LineStrip new_linestrip(*vertex, linestrip_it.iterator->template get_end<reverse>());
    LineStrip new_linestrip(vertex.parent->get_parent());
    new_linestrip.set_start<reverse>(*vertex);
    new_linestrip.set_end<reverse>(linestrip_it.iterator->template get_end<reverse>());

    linestrip_it.iterator->template set_end<reverse>(*vertex);
    linestrip_it.insert_after_and_move_to_new(new_linestrip, parent->linestrips);

    vertex.parent = &(*linestrip_it.iterator);
    vertex.set_to_start();

    /* old */
    /*
    linestrip_iterator<reverse> insert_position = linestrip_it;
    ++ insert_position.iterator;
    parent->linestrips.insert(insert_position.base(), new_linestrip);
    -- insert_position.iterator;

    linestrip_it.iterator = insert_position.iterator;
    vertex.set_to_start_of(&*linestrip_it.iterator);
    */
}


templates
void LSS::Vertex<reverse>::set_linestrip_start(float alpha)
{
    EdgePoint new_start = to_edge_point(alpha);
    linestrip_it.iterator->template set_start<reverse>(new_start);
}

templates
bool LSS::Vertex<reverse>::at_end()
{
    return (linestrip_it.at_end_of(parent->linestrips) && vertex.at_end());
}
templates
bool LSS::Vertex<reverse>::at_start()
{
    return (linestrip_it.at_start_of(parent->linestrips) && vertex.at_start());
}
templates
void LSS::Vertex<reverse>::set_to_start()
{
    linestrip_it.set_to_start_of(parent->linestrips);
    vertex = LineStrip::Vertex<reverse>(&(*linestrip_it.iterator));
    vertex.set_to_start();
}

templates
LSS::Vertex<reverse>& LSS::Vertex<reverse>::operator ++()
{
    ++ vertex;
    if (vertex.at_end()) {
        if ( ! linestrip_it.at_end_of(parent->linestrips)) {
            ++ linestrip_it.iterator; 
            vertex.set_to_start_of(&*linestrip_it.iterator);
        }
    }
    return *this;
}

templates
LSS::Vertex<reverse>& LSS::Vertex<reverse>::operator --()
{
    if (vertex.at_start()) {
        if ( ! linestrip_it.at_start_of(parent->linestrips)) {
            -- linestrip_it.iterator; 

            vertex = LineStrip::Vertex<reverse>(&(*linestrip_it.iterator));
            vertex.set_to_end();
        }
    } else {
        -- vertex;
    }
    return *this;
}


templates
glm::vec2 LSS::Vertex<reverse>::entry_point()
{
    if (vertex.at_start() && ! linestrip_it.at_start_of(parent->linestrips)) {
        linestrip_iterator<reverse> previous_linestrip = linestrip_it;
        -- previous_linestrip.iterator;
        return parent->align_transform * previous_linestrip.iterator->template get_end<reverse>().point_t();
    } else {
        return parent->align_transform * vertex.point_t();
    }
}
templates
glm::vec2 LSS::Vertex<reverse>::exit_point()
{
    return parent->align_transform * vertex.point_t();
}

templates
glm::vec2 LSS::Vertex<reverse>::best_point()
{
    glm::vec2 entry = entry_point();
    glm::vec2 exit = exit_point();
    if (x_bias == LEFT) {
        if (entry.x < exit.x)
            return entry;
        else
            return exit;
        
    } else { //RIGHT
        if (entry.x > exit.x)
            return entry;
        else
            return exit;
    }
}

#undef LSS
#undef templates
