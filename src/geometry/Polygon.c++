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
#include <cmath>
#include <stdexcept>
#include <set>
#include <sstream>

#define DRAW_VERTEX_NUMBERS true

/*
Some useful sources I used. (not exhaustive)
http://www.cs.uu.nl/docs/vakken/ga/slides3.pdf
http://www.cs.unc.edu/~dm/CODE/GEM/chapter.html#Fournier84
https://www.cs.ucsb.edu/~suri/cs235/Triangulation.pdf

*/
extern FontRenderer *fontRenderer;

int abs_mod(int n, int range)
{
    return (n + range) % range;
}

Polygon::Polygon()
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
float Polygon::radius()
{
	glm::vec2 c = centroid();
	// Find point whose distance is MAXIMUM from c.
	float m = 0;
	float d;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		d = distance(c, *it); // TODO: made this relative to local coor system
		if (d > m) m = d;
	}
	return m;
}

glm::vec2 Polygon::transform(glm::vec2 point)
{
    return glm::vec2(   matrix[0][0] * point.x + matrix[1][0] * point.y + matrix[2][0],
                        matrix[0][1] * point.x + matrix[1][1] * point.y + matrix[2][1] );
}
glm::vec2 Polygon::transformed(int vertex_index)
{
    glm::vec2 point = vertices[vertex_index];
    return glm::vec2(   matrix[0][0] * point.x + matrix[1][0] * point.y + matrix[2][0],
                        matrix[0][1] * point.x + matrix[1][1] * point.y + matrix[2][1] );
}
/* translates first such that the center is at origin */
glm::vec2 Polygon::transform_center(glm::vec2 point, glm::vec2 center)
{
    point -= center;
    return glm::vec2(   matrix[0][0] * point.x + matrix[1][0] * point.y + matrix[2][0],
                        matrix[0][1] * point.x + matrix[1][1] * point.y + matrix[2][1] );
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

// TRIANGLE FAN TODO make this clear
//void Polygon::append_stencil_triangles(std::vector<float> &buffer)
void Polygon::append_stencil_triangles(BufferWriter<float> &buffer)
{
    // i is the edge index
    for (uint i = 0; i < vertices.size(); i ++)
    {
        buffer.write(vertices[i].x, vertices[i].y);
        glm::vec2 transformed = transform(vertices[i]);
        if (DRAW_VERTEX_NUMBERS) {
            fontRenderer->setColor(255, 255, 255);
            fontRenderer->addText(std::to_string(i), transformed.x, transformed.y,  false);
        }
    }
}
void Polygon::append_lines_to_vector(std::vector<float> &list)
{
    for (uint i = 0; i < vertices.size(); i ++) {
        int j = i + 1; j %= vertices.size();
        glm::vec2 vec_i = transform(vertices[i]);
        glm::vec2 vec_j = transform(vertices[j]);
        list.push_back(vec_i.x);
        list.push_back(vec_i.y);
        list.push_back(vec_j.x);
        list.push_back(vec_j.y);
    }
}

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
    return parent->transform(parent->vertices[index]);
}

//////////////////////////////////////////////
//
//              Edge pointer
//
//////////////////////////////////////////////

Polygon::Edge::Edge(int index, Polygon *parent)
    :index(index % parent->vertices.size()), parent(parent) { }

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

bool Polygon::Edge::operator== (Edge other)
{
    return index == other.index && parent == other.parent;
}
