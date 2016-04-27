#include "LineStrip.h"
void LineStrip::append_lines_to_vector(std::vector<float> &list, float r, float g, float b)
{
    Vertex<false> it(Vertex<false>::START_INDEX, this);
    Vertex<false> prev = it;
    while (true) {
        glm::vec2 vec_i = prev.point_t();
        glm::vec2 vec_j = it.point_t();
        list.push_back(vec_i.x);
        list.push_back(vec_i.y);
        list.push_back(r); list.push_back(g); list.push_back(b);
        list.push_back(vec_j.x);
        list.push_back(vec_j.y);
        list.push_back(r); list.push_back(g); list.push_back(b);
        if (it.at_end()) break;
        prev = it;
        ++ it;
    }
}

const bool LineStrip::operator== (const LineStrip& other)
{
    return (start == other.start && end == other.end && parent == other.parent);
}

///////////////////
/*** LineStrip ***/
///////////////////

template <>
void LineStrip::set_start<false>(EdgePoint value) { start = value; }

template <>
void LineStrip::set_end<false>(EdgePoint value) { end = value; }

template <>
void LineStrip::set_start<true>(EdgePoint value) { end = value; }

template <>
void LineStrip::set_end<true>(EdgePoint value) { start = value; }

template<>
EdgePoint LineStrip::get_start<false>() { return start; }

template<>
EdgePoint LineStrip::get_end<false>() { return end; }

template<>
EdgePoint LineStrip::get_start<true>() { return end; }

template<>
EdgePoint LineStrip::get_end<true>() { return start; }

///////////////////////////
/*** LineStrip::Vertex ***/
///////////////////////////
template <>
EdgePoint LineStrip::Vertex<false>::to_edge_point(float alpha)
{
    if (index == START_INDEX) {
        if (parent->start.index == parent->end.index) {
            // START to END
            return EdgePoint(parent->start.index, parent->start.alpha + (parent->end.alpha - parent->start.alpha) * alpha, parent->parent);
        } else {
            // START to INTERMEDIATE
            return EdgePoint(parent->start.index, parent->start.alpha + alpha - parent->start.alpha * alpha, parent->parent);
        }
    } else if (index == END_INDEX) {
        if (alpha == 0) {
            return parent->end;
        } else {
            assert(!"Trying to move past END_INDEX in LineStrip::Vertex::to_edge_point");
        }
    } else {
        if (index == parent->end.index) { 
            // INTERMEDIATE to END
            return EdgePoint(index, parent->end.alpha * alpha, parent->parent);
        } else {
            // INTERMEDIATE to INTERMEDIATE
            return EdgePoint(index, alpha, parent->parent);
        }
    }
}

// where alpha represents progress in reverse direction from this LineStrip::Vertex
template <>
EdgePoint LineStrip::Vertex<true>::to_edge_point(float alpha)
{
    if (index == START_INDEX) {
        if (alpha == 0) {
            return parent->start;
        } else {
            assert(!"Trying to move past START_INDEX in LineStrip::Vertex::to_edge_point_reverse");
        }
    } else if (index == END_INDEX) {
        if (parent->end.index == parent->start.index) {
            // END to START
            return EdgePoint(parent->end.index,
                    parent->end.alpha - (parent->end.alpha - parent->start.alpha) * alpha,
                    parent->parent);
        } else {
            // END to INTERMEDIATE
            return EdgePoint(parent->end.index,
                    parent->end.alpha - parent->end.alpha * alpha,
                    parent->parent);
        }
    } else {
        int next_index = index - 1;
        if (index < 0) index += parent->parent->vertices.size();
        if (next_index == parent->start.index) {
            // INTERMEDIATE to START
            return EdgePoint(index - 1, parent->start.alpha + (1 - parent->start.alpha) * (1 - alpha), parent->parent);
        } else {
            // INTERMEDIATE to INTERMEDIATE
            return EdgePoint(index - 1, 1 - alpha, parent->parent);
        }
    }
}

std::ostream& operator<< (std::ostream& out, LineStrip& ls)
{
    out << "(" << ls.get_start<false>().index << ", " << ls.get_start<false>().alpha << ") -> (" <<
        ls.get_end<false>().index << ", " << ls.get_end<false>().alpha << ")";
    return out;
}
