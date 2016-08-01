#include "Polygon.h"
#include "geometry.h"

#include "tmp.h"
#include "LinkedList.h"
#include "glutils.h"
#include "render/font/FontRenderer.h"

#include <list>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <set>
#include <sstream>
#include <limits>

#define DRAW_VERTEX_NUMBERS false

extern FontRenderer *fontRenderer;

int abs_mod(int n, int range)
{
    return (n + range) % range;
}

Polygon::Polygon()
    : orientation{}, moment_of_inertia{}, mass{}, center_of_mass{}, CCW{}
{
}

Polygon::Polygon(std::vector<glm::vec2> vertices)
{
    _vertices = vertices;
    calculate_shape_dependent_variables();
}


float Polygon::signed_area()
{
	float A = 0;
	std::vector<glm::vec2>::iterator next;
	for (auto it = _vertices.begin(); it != _vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == _vertices.end())
            next = _vertices.begin();
		A += it->x * next->y - next->x * it->y;
	}
	A *= 0.5f;
	return A;
}

void Polygon::calculate_shape_dependent_variables()
{
	float A = signed_area();
	glm::vec2 C {};
    glm::vec2 a, b;
	std::vector<glm::vec2>::iterator next;
	for (auto it = _vertices.begin(); it != _vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == _vertices.end())
            next = _vertices.begin();
		C.x += (it->x + next->x) * (it->x * next->y - next->x * it->y);
		C.y += (it->y + next->y) * (it->x * next->y - next->x * it->y);
	}
	C /= 6.f * A;

    this->mass = fabs(A); 
    this->moment_of_inertia = calculate_moment_of_inertia();
    this->CCW = (A > 0);
    this->center_of_mass = C;
    apply_center_of_mass(center_of_mass);
    this->radius = calculate_radius();
}
float Polygon::calculate_moment_of_inertia()
{
    float sum1=0;
    float sum2=0;
    std::vector<glm::vec2>::iterator next;
	for (auto it = _vertices.begin(); it != _vertices.end(); it ++)
	{
		next = it; next ++;
        if (next == _vertices.end())
            next = _vertices.begin();
        sum1 += cross(*next, *it) * (glm::dot(*next, *next) + glm::dot(*next, *it) + glm::dot(*it, *it));
        sum2 += cross(*next, *it);
    }
    return (mass/6*sum1/sum2);
}
float Polygon::calculate_radius()
{
	// Find point whose distance is MAXIMUM from origin - assuming we have already centered the center of mass
	float m = 0;
	for (auto it = _vertices.begin(); it != _vertices.end(); it ++)
	{
		float d = length_squared(*it);
		if (d > m) m = d;
	}
	return sqrt(m);
}
AABB<2> Polygon::calc_bounding_box()
{
    AABB<2> box(std::numeric_limits<float>::max(), - std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max(), - std::numeric_limits<float>::max());
    //TODO set to MIN and MAX
	for (::Vertex v : vertices())
    {
        box.min[0] = std::min(box.min[0], v.point.x);
        box.min[1] = std::min(box.min[1], v.point.y);
        box.max[0] = std::max(box.max[0], v.point.x);
        box.max[1] = std::max(box.max[1], v.point.y);
    }
    return box;
}

bool Polygon::vertex_is_concave(int index)
{
    glm::vec2& prev = _vertices[(index - 1 + _vertices.size())   % _vertices.size()];
    glm::vec2& next = _vertices[(index + 1)                      % _vertices.size()];
    glm::vec2& vert = _vertices[index];
    return leftof(next - prev, vert - prev) ^ CCW;
}

void Polygon::apply_center_of_mass(glm::vec2 center)
{
    for (int i = 0; i < _vertices.size(); i ++) {
        _vertices[i] -= center;
    }
}

glm::vec2 Polygon::transform(glm::vec2 point)
{
    // return glm::vec2(   matrix[0][0] * point.x + matrix[1][0] * point.y + matrix[2][0],
                        // matrix[0][1] * point.x + matrix[1][1] * point.y + matrix[2][1] );
    float c = cos(orientation);
    float s = sin(orientation);
    return glm::vec2(   c * point.x - s * point.y + position.x,
                        s * point.x + c * point.y + position.y);
}
glm::vec2 Polygon::detransform(glm::vec2 point)
{
    // translate back, rotate back
    point -= position;
    float c = cos(orientation);
    float s = sin(orientation);
    return glm::vec2(   + c * point.x + s * point.y,
                        - s * point.x + c * point.y);
}
glm::vec2 Polygon::transformed(int vertex_index)
{
    vertex_index = (vertex_index + _vertices.size()) % _vertices.size();
    glm::vec2 point = _vertices[vertex_index];
    return transform(point);
}
/* translates first such that the center is at origin */
glm::vec2 Polygon::transform_center(glm::vec2 point, glm::vec2 center)
{
    point -= center;
    return transform(point);
}

glm::vec2 Polygon::edge_point(int index, float alpha)
{
    int new_index = (index + 1) % _vertices.size();
    return (1 - alpha) * transformed(index) + alpha * transformed(new_index);
}
glm::vec2 Polygon::model_edge_point(int index, float alpha)
{
    int new_index = (index + 1) % _vertices.size();
    return (1 - alpha) * _vertices[index] + alpha * _vertices[new_index];
}

/* 
AABB<2> Polygon::calculate_bounding_box()
{
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::min();
    float max_y = std::numeric_limits<float>::min();
    for (int i = 0; i < _vertices.size(); i ++) {
        glm::vec2 point = transformed(i);
        if (point.x < min_x) min_x = point.x;
        if (point.x < max_x) max_x = point.x;
        if (point.y < min_y) min_y = point.y;
        if (point.y < max_y) max_y = point.y;
    }
    return AABB<2>(min_x, max_x, min_y, max_y);
} */


//TODO WARNING: using a set for intersecting edges, edges with same y-value of start vertex is not allowed!

// Sorting: could sort a list of indices but complicates the comparision function

std::ostream &operator << (std::ostream &lhs, glm::vec2 &rhs);


int Polygon::num_edges()
{
	return _vertices.size();
}
//TODO Caution: Local coordinates
LineSegment Polygon::get_edge(int index)
{
	if (index < 0 || index >= num_edges()) throw(std::out_of_range("Edge out of range."));
	glm::vec2 s = _vertices[index];
	glm::vec2 t = _vertices[(index < num_edges() - 1)? index + 1 : 0];
	// s += position;
	// t += position;
	return LineSegment(s, t);
}

//////////////////////////////////////////////
//
//              Vertex pointer
//
//////////////////////////////////////////////

Polygon::Vertex::Vertex(int index, Polygon* parent)
    :index(index % parent->_vertices.size()), parent(parent) { }

void Polygon::Vertex::set_index(int val)
{
    index = val;
    if (index < 0) index += parent->_vertices.size();
    index %= parent->_vertices.size();
}
glm::vec2& Polygon::Vertex::operator* ()
{
	return parent->_vertices[index];
}
glm::vec2* Polygon::Vertex::operator-> ()
{
	return &(parent->_vertices[index]);
}
glm::vec2& Polygon::Vertex::successive()
{
	int i = (index == parent->_vertices.size() - 1) ? 0 : index + 1;
	return parent->_vertices[i];
}
glm::vec2& Polygon::Vertex::preceding()
{
	int i = (index == 0) ? parent->_vertices.size() - 1 : index - 1;
	return *(&(parent->_vertices[i]));
}
Polygon::Vertex& Polygon::Vertex::operator++ ()
{
	// index = (index == parent->_vertices.size() - 1) ? 0 : index + 1;
    index ++;
    index %= parent->_vertices.size();
    return *this;
}
Polygon::Vertex& Polygon::Vertex::operator-- ()
{
	// index = (index == 0) ? parent->_vertices.size() - 1 : index - 1;
    index --;
    if (index < 0) index += parent->_vertices.size();
    return *this;
}
bool Polygon::Vertex::operator== (Polygon::Vertex& v)
{
    return (index == v.index && parent == v.parent);
}
glm::vec2 Polygon::Vertex::transformed()
{
    glm::vec2 result = parent->transformed(index);
    return result;
}

//////////////////////////////////////////////
//
//              Edge pointer
//
//////////////////////////////////////////////

Polygon::Edge::Edge(int index, Polygon *parent)
    : parent(parent)
{
    this->index = (index + parent->_vertices.size()) % parent->_vertices.size();
}

Polygon::Edge Polygon::first_edge() {
    return Edge(0, this);
}
Polygon::Edge Polygon::last_edge() {
    return -- first_edge();
}
glm::vec2 Polygon::Edge::normal_tr()
{
    // Needs to know CCW of polygon to know the correct direction (+/-)
    glm::vec2 edge_vec = end_tr() - start_tr();
    if (parent->CCW) {
        /* Normal always points out of the polygon */
        edge_vec = - edge_vec;
    }
    return glm::normalize(glm::vec2( - edge_vec.y, edge_vec.x));
}

glm::vec2& Polygon::Edge::start() const
{
	return parent->_vertices[index];
}
glm::vec2& Polygon::Edge::end() const
{
	int i = (index == parent->_vertices.size() - 1) ? 0 : index + 1;
	return parent->_vertices[i];
}
glm::vec2 Polygon::Edge::start_tr() const
{
	return parent->transformed(index);
}
glm::vec2 Polygon::Edge::end_tr() const
{
	int i = (index == parent->_vertices.size() - 1) ? 0 : index + 1;
	return parent->transformed(i);
}

// Find y value given x value. Doesn't care about bounds.
int Polygon::Edge::operator() (int x) const
{
	glm::vec2 delta = end() - start();
	assert(start().x != end().x);
	float t = (x - start().x)/delta.x;
	return start().y + t * delta.y;
}
bool Polygon::Edge::operator== (Polygon::Edge other) {
    return index == other.index;
}
bool Polygon::Edge::operator!= (Polygon::Edge other) {
    return index != other.index;
}
// Experimental iterator...
Polygon::Edge& Polygon::Edge::operator++ ()
{
	 index = (index == parent->_vertices.size() - 1) ? 0 : index + 1;
     return *this;
}
Polygon::Edge& Polygon::Edge::operator-- ()
{
    index = (index == 0) ? (parent->_vertices.size() - 1) : index - 1;
    return *this;
}


