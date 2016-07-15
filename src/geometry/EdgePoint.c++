#include "EdgePoint.h"
#include "Polygon.h"
#include "geometry.h"

EdgePoint::EdgePoint(int index, float alpha, Polygon* parent)
    :alpha(alpha), parent(parent)
{
    if (index < 0) index += parent->num_vertices();
    this->index = index % parent->num_vertices();
}
glm::vec2 EdgePoint::point()
{
    int next_index = (index + 1) % parent->num_vertices();
    return (1 - alpha) * parent->model_vertex(index) + alpha * parent->model_vertex(next_index);
}
glm::vec2 EdgePoint::point_t()
{
    int next_index = (index + 1) % parent->num_vertices();
    return parent->transform(parent->model_vertex(index) + alpha * (parent->model_vertex(next_index) - parent->model_vertex(index)));
}
std::ostream& operator<< (std::ostream& o, EdgePoint p){
    o << "(" << p.index << ", " << p.alpha << ", " << p.parent << ")";
    return o;
}
const bool EdgePoint::operator== (const EdgePoint& other)
{
    return index == other.index && alpha == other.alpha && parent == other.parent;
}
