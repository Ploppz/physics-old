#include "Polygon.h"
#include "geometry.h"

#include "../tmp.h"
#include "../LinkedList.h"
#include "../glutils.h"
#include "../typewriter/FontRenderer.h"

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

#define DRAW_VERTEX_NUMBERS false

extern FontRenderer *fontRenderer;

int abs_mod(int n, int range)
{
    return (n + range) % range;
}

Polygon::Polygon()
    : orientation{}, moment_of_inertia{}, mass{}, center_of_mass{}, CCW{} // Why is this -nan if not initialized??
{
}


float Polygon::signed_area()
{
	float A = 0;
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();
		A += it->x * next->y - next->x * it->y;
	}
	A *= 0.5f;
	return A;
}

// LOCAL
glm::vec2 Polygon::centroid()
{
	float A = signed_area();
	glm::vec2 C {};
    glm::vec2 a, b;
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();
		C.x += (it->x + next->x) * (it->x * next->y - next->x * it->y);
		C.y += (it->y + next->y) * (it->x * next->y - next->x * it->y);
	}
	C /= 6.f * A;

	return C;
}
void Polygon::calculate_shape_dependent_variables()
{
	float A = signed_area();
	glm::vec2 C {};
    glm::vec2 a, b;
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();
		C.x += (it->x + next->x) * (it->x * next->y - next->x * it->y);
		C.y += (it->y + next->y) * (it->x * next->y - next->x * it->y);
	}
	C /= 6.f * A;

    mass = fabs(A); 
    CCW = (A > 0);
    center_of_mass = C;
    std::cout << "CCW: " << CCW << std::endl;

    moment_of_inertia = calculate_moment_of_inertia();
    std::cout << "MOMENT : " << moment_of_inertia << std::endl; 
}
float lengthcross (glm::vec2 a, glm::vec2 b) {
   return (fabs(a.x*b.y - a.y*b.x));
}
float Polygon::calculate_moment_of_inertia()
{
    float sum1=0;
    float sum2=0;
    std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
        if (next == vertices.end())
            next = vertices.begin();
        sum1 += lengthcross(*next, *it) * (glm::dot(*next, *next) + glm::dot(*next, *it) + glm::dot(*it, *it));
        sum2 += lengthcross(*next, *it);
    }
    return (mass/6*sum1/sum2);
}
float Polygon::radius()
{
	glm::vec2 c = centroid();
	// Find point whose distance is MAXIMUM from c.
	float m = 0;
	float d;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		d = distance(c, *it);
		if (d > m) m = d;
	}
	return m;
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
glm::vec2 Polygon::transformed(int vertex_index)
{
    vertex_index = (vertex_index + vertices.size()) % vertices.size();
    glm::vec2 point = vertices[vertex_index];
    return transform(point);
}
/* translates first such that the center is at origin */
glm::vec2 Polygon::transform_center(glm::vec2 point, glm::vec2 center)
{
    point -= center;
    return transform(point);
}
glm::vec2 Polygon::get_point(int vertex_number, float alpha)
{
    int next_vertex_number = vertex_number + 1;
    if (next_vertex_number >= vertices.size()) next_vertex_number = 0;
    return  transform(vertices[vertex_number]) * (1 - alpha) + transform(vertices[next_vertex_number]) * alpha;
}


//TODO WARNING: using a set for intersecting edges, edges with same y-value of start vertex is not allowed!

// Sorting: could sort a list of indices but complicates the comparision function

std::ostream &operator << (std::ostream &lhs, glm::vec2 &rhs);


int Polygon::num_edges()
{
	return vertices.size();
}
//TODO Caution: Local coordinates
LineSegment Polygon::get_edge(int index)
{
	if (index < 0 || index >= num_edges()) throw(std::out_of_range("Edge out of range."));
	glm::vec2 s = vertices[index];
	glm::vec2 t = vertices[(index < num_edges() - 1)? index + 1 : 0];
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
    :index(index % parent->vertices.size()), parent(parent) { }

void Polygon::Vertex::set_index(int val)
{
    index = val;
    if (index < 0) index += parent->vertices.size();
    index %= parent->vertices.size();
}
glm::vec2& Polygon::Vertex::operator* ()
{
	return parent->vertices[index];
}
glm::vec2* Polygon::Vertex::operator-> ()
{
	return &(parent->vertices[index]);
}
glm::vec2& Polygon::Vertex::successive()
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}
glm::vec2& Polygon::Vertex::preceding()
{
	int i = (index == 0) ? parent->vertices.size() - 1 : index - 1;
	return *(&(parent->vertices[i]));
}
Polygon::Vertex& Polygon::Vertex::operator++ ()
{
	// index = (index == parent->vertices.size() - 1) ? 0 : index + 1;
    index ++;
    index %= parent->vertices.size();
    return *this;
}
Polygon::Vertex& Polygon::Vertex::operator-- ()
{
	// index = (index == 0) ? parent->vertices.size() - 1 : index - 1;
    index --;
    if (index < 0) index += parent->vertices.size();
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
    :index(index % parent->vertices.size()), parent(parent) { }

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
	return parent->vertices[index];
}
glm::vec2& Polygon::Edge::end() const
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}
glm::vec2 Polygon::Edge::start_tr() const
{
	return parent->transform(parent->vertices[index]);
}
glm::vec2 Polygon::Edge::end_tr() const
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->transform(parent->vertices[i]);
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
	 index = (index == parent->vertices.size() - 1) ? 0 : index + 1;
     return *this;
}
Polygon::Edge& Polygon::Edge::operator-- ()
{
    index = (index == 0) ? (parent->vertices.size() - 1) : index - 1;
    return *this;
}
