#include "EdgePoint.h"
#include "Polygon.h"

EdgePoint::EdgePoint(int index, float alpha, Polygon* parent)
    :alpha(alpha), parent(parent)
{
    if (index < 0) index += parent->vertices.size();
    this->index = index % parent->vertices.size();
}
glm::vec2 EdgePoint::point()
{
    int next_index = (index + 1) % parent->vertices.size();
    return (1 - alpha) * parent->vertices[index] + alpha * parent->vertices[next_index];
}
glm::vec2 EdgePoint::point_t()
{
    int next_index = (index + 1) % parent->vertices.size();
    // return parent->transform((1 - alpha) * parent->vertices[index] + alpha * parent->vertices[next_index]);
    return parent->transform(parent->vertices[index] + alpha * (parent->vertices[next_index] - parent->vertices[index]));
}
std::ostream& operator<< (std::ostream& o, EdgePoint p){
    o << "(" << p.index << ", " << p.alpha << ", " << p.parent << ")";
    return o;
}
const bool EdgePoint::operator== (const EdgePoint& other)
{
    return index == other.index && alpha == other.alpha && parent == other.parent;
}
