#pragma once
#include <iostream>
#include <vector>
#include "../constants.h"
#include "../geometry.h"
#include "../Intersection.h"
#include <glm/glm.hpp>
/* External classes */
class Polygon;
/**/

/**
Revers or not:
Internal representation is the same, it's just different
ways to work with the LineStrip.
-
if reverse:
    the interface virtually swaps start and end.
    Vertex iterates backwards.
**/

class LineStrip
{
 public: // Classes
    template <bool reverse>
    class Vertex;
 public:
    LineStrip(Polygon* parent) : start(0, 0, parent), end(0, 0, parent), parent(parent)
    {}
    LineStrip(EdgePoint start, EdgePoint end) : start(start), end(end), parent(start.parent)
    {
        assert(start.parent == end.parent);
    } 

    void swap_start_and_end();

    template <bool reverse>
    void set_start(const EdgePoint value);

    template <bool reverse>
    void set_end(const EdgePoint value);

    template <bool reverse>
    EdgePoint get_start();

    template <bool reverse>
    EdgePoint get_end();

    Polygon* get_parent() {return parent;}
    
    // bool CCW;

    void append_lines_to_vector(std::vector<float> &list, float r, float g, float b);

    const bool operator== (const LineStrip& other);

 private:
    EdgePoint start, end;
    Polygon* parent; // each EdgePoint also has one - questionable if needed
};

std::ostream& operator<< (std::ostream& out, LineStrip& ls);

template <bool reverse>
class LineStrip::Vertex {
 private:
    int index;
    void increment();
    void decrement();
 public:
    Vertex(): index(START_INDEX), parent(0) {}; 
    Vertex(LineStrip *parent);
    Vertex(int index, LineStrip *parent);
    int get_index() { return index; }
    float get_alpha();
    LineStrip* get_parent() { return parent; }
    bool at_start();
    bool at_end();
    void set_to_start();
    void set_to_end();
    void set_index(int val);
    void set_to_start_of(LineStrip* new_parent);
    // A point between two edge points is a new edge point...
    EdgePoint to_edge_point(float alpha);
    
    EdgePoint operator* ();

    // Changes what vertex is pointed to
    Vertex& operator++ ();
    Vertex& operator-- ();
    //
    bool operator== (Vertex& v);
    glm::vec2 point();
    glm::vec2 point_t(); //transformed point

    LineStrip *parent;
public:
    static const int START_INDEX = -1;
    static const int END_INDEX = -2;
};

/** Index (helper) **/
/** Helper struct for templaping reverse logic. **/
template <bool reverse>
struct Index;
template<>
struct Index<false> {
    static const int START = LineStrip::Vertex<false>::START_INDEX;
    static const int END = LineStrip::Vertex<false>::END_INDEX;
};
template<>
struct Index<true> {
    static const int START = LineStrip::Vertex<false>::END_INDEX;
    static const int END = LineStrip::Vertex<false>::START_INDEX;
};
///////////////////////////
/*** LineStrip::Vertex ***/
///////////////////////////

// TODO dunno if I should specialize for reverse=true
template <bool reverse>
LineStrip::Vertex<reverse>::Vertex(int index, LineStrip* parent)
    :index(index % (int)parent->parent->vertices.size()), parent(parent)
{
}

template <bool reverse>
LineStrip::Vertex<reverse>::Vertex(LineStrip *parent)
    :index(START_INDEX), parent(parent)
{
}


template <bool reverse>
bool LineStrip::Vertex<reverse>::at_start()
{
    return index == Index<reverse>::START;
}

template <bool reverse>
bool LineStrip::Vertex<reverse>::at_end()
{
    return index == Index<reverse>::END;
}

template <bool reverse>
void LineStrip::Vertex<reverse>::set_to_start()
{
    index = Index<reverse>::START;
}

template <bool reverse>
void LineStrip::Vertex<reverse>::set_to_end()
{
    index = Index<reverse>::END;
}

template <bool reverse>
float LineStrip::Vertex<reverse>::get_alpha() {
    if (index == START_INDEX) return parent->start.alpha;
    else if (index == END_INDEX) return parent->end.alpha;
    else return 0;
}

template <bool reverse>
void LineStrip::Vertex<reverse>::set_to_start_of(LineStrip* new_parent)
{
    parent = new_parent;
    set_to_start();
}
template <bool reverse>
glm::vec2 LineStrip::Vertex<reverse>::point()
{
    if (index == START_INDEX) {
        return parent->start.point();
    } else if (index == END_INDEX) {
        return parent->end.point();
    } else {
        return parent->parent->vertices[index];
    }
}
template <bool reverse>
glm::vec2 LineStrip::Vertex<reverse>::point_t()
{
    if (index == START_INDEX) {
        return parent->parent->transform( parent->start.point() );
    } else if (index == END_INDEX) {
        return parent->parent->transform( parent->end.point() );
    } else {
        return parent->parent->transform( parent->parent->vertices[index] );
    }
}

template <bool reverse>
EdgePoint LineStrip::Vertex<reverse>::operator* ()
{
    if (index == START_INDEX) {
        return parent->start;
    } else if (index == END_INDEX) {
        return parent->end;
    } else {
        return EdgePoint(index, 0, parent->parent);
    }
}
template <bool reverse>
void LineStrip::Vertex<reverse>::increment()
{
    if (index == START_INDEX) {
        if (parent->start.index == parent->end.index)
            index = END_INDEX;       
        else
            index = parent->start.index; // Note to self: not +1, because that comes in the next if!
    }
    if (index != END_INDEX) {
        if (index == parent->end.index) {
            index = END_INDEX;
        } else {
            index ++;
            index %= parent->parent->vertices.size();
            if (index == parent->end.index && parent->end.alpha == 0) {
                // In this case, the next index was indeed the end
                index = END_INDEX;
            }
        }
    }
}

template <bool reverse>
void LineStrip::Vertex<reverse>::decrement ()
{
    if (index == END_INDEX) {
        if (parent->end.index == parent->start.index)
            index = START_INDEX;
        else
            if (parent->end.alpha == 0) {
                index = parent->end.index - 1;
                if (index < 0)  index += parent->parent->vertices.size();
                if (index == parent->start.index)
                    index = START_INDEX;
            } else {
                index = parent->end.index;
            }
    }
    else if (index != START_INDEX) {
        index --;
        if (index < 0)
            index += parent->parent->vertices.size();
        if (index == parent->start.index)
            index = START_INDEX;
    }
}
template <bool reverse>
typename LineStrip::Vertex<reverse>& LineStrip::Vertex<reverse>::operator++ ()
{
    if (reverse) {
        decrement();
    } else {
        increment();
    }
    return *this;
}
template <bool reverse>
typename LineStrip::Vertex<reverse>& LineStrip::Vertex<reverse>::operator-- ()
{
    if (reverse) {
        increment();
    } else {
        decrement();
    }
    return *this;
}
template <bool reverse>
bool LineStrip::Vertex<reverse>::operator== (LineStrip::Vertex<reverse>& v)
{
    return (index == v.index && parent == v.parent);
}
