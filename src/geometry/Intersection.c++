#include "Polygon.h"
#include "Intersection.h"
#include "geometry.h"
#include "../LinkedList.h"
#include "../Renderer.h"
#include "../tmp.h"
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <limits>
#include <utility>

using namespace glm;



////////////////////
// Intersection
///////////////////
Intersect::Intersect(Polygon::Edge edge1, Polygon::Edge edge2)
{
    this->edge1 = edge1; this->edge2 = edge2;
    // Find intersection point (last argument is output)
    bool result = intersect(edge1.start_tr(), edge1.end_tr(),
            edge2.start_tr(), edge2.end_tr(), point, alpha1, alpha2);
    if (!result) {
        std::cout << "Error: edges do not actually overlap (" << edge1.get_index() << ", " << edge2.get_index() << ")" << std::endl;
    }
}




/****************/
/* HybridVertex */
/****************/
HybridVertex::HybridVertex(Intersect& i)
{
    edge1_owner = i.edge1.get_parent();
    edge1_index = i.edge1.get_index();
    edge2_owner = i.edge2.get_parent();
    edge2_index = i.edge2.get_index();
    intersect = true;
    //
    point = i.point;
}
HybridVertex::HybridVertex(Polygon::Vertex v)
{
    owner = v.get_parent();
    vertex = v.get_index();
    intersect = false;
    point = v.transformed();
}
/****************/
/* Intersection */
/****************/
Intersection::Intersection(Polygon& p)
{
    vertices.resize(p.vertices.size());
    for (int i = 0; i < p.vertices.size(); i ++)
    {
        HybridVertex v( Polygon::Vertex(i, &p) );
        vertices[i] = v;
    }
}

float Intersection::signed_area()
{
	float A = 0;
	std::vector<HybridVertex>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();
		A += it->point.x * next->point.y - next->point.x * it->point.y;
	}
	A *= 0.5f;
	return A;
}
glm::vec2 Intersection::centroid()
{
	float A = signed_area();
	glm::vec2 C {};
    glm::vec2 a, b;
	std::vector<HybridVertex>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();
		C.x += (it->point.x + next->point.x) * (it->point.x * next->point.y - next->point.x * it->point.y);
		C.y += (it->point.y + next->point.y) * (it->point.x * next->point.y - next->point.x * it->point.y);
	}
	C /= 6.f * A;

	return C;
}

glm::vec2 Intersection::find_normal_towards(Polygon* p) // or.. find direction of intersect wrt
{
    // PERFORMANCE does a lot of similar things as find_depth etc,
    int intersect1 = -1, intersect2;
    for (int i = 0; i < vertices.size(); i ++) {
        if (vertices[i].intersect) {
            if (intersect1 == -1) {
                intersect1 = i;
            } else {
                intersect2 = i;
            }
        }
    }
    /* PERFORMANCE don't know whether we need to normalize */
    glm::vec2 normal = glm::normalize(vertices[intersect2].point - vertices[intersect1].point);
    normal = glm::vec2( - normal.y, normal.x );
    // normal points to the left 
    // found with assumption that Polygon p lies to the left of (intersect2 - intersect1)
    /* check if polygon is left of the line .. */
    glm::vec2 next_edge = Vertex(intersect1, this).successive().point - Vertex(intersect1, this)->point;
    glm::vec2 prev_edge = Vertex(intersect1, this).preceding().point - Vertex(intersect1, this)->point;
    if (leftof(next_edge, prev_edge)) {
        if (edge_owner(intersect1) != p) {
            normal = - normal;
        }
    } else {
        if (edge_owner(intersect1) == p) {
            normal = - normal;
        }
    }
    return normal;
}
Polygon* Intersection::edge_owner(int edge_start_index)
{
    HybridVertex v1 = vertices[edge_start_index];
    HybridVertex v2 = vertices[clamp_index(edge_start_index + 1)];

    if ( ! v1.intersect) {
        return v1.owner;
    } else if ( ! v2.intersect) {
        return v2.owner;
    } else {
        /* Find the polygon whose _same single edge_ is present in the two intersections */
        if (v1.edge1_index == v2.edge1_index || v1.edge1_index == v2.edge2_index) {
            return v1.edge1_owner;
        } else {
            /* Fair assumption that now, v1.edge2_index must be the same as any of v2's edges */
            assert (v1.edge2_index == v2.edge1_index || v1.edge2_index == v2.edge2_index);
            return v1.edge2_owner;
        }
    }
}
int Intersection::clamp_index(int index)
{
    if (index < 0) index += vertices.size();
    index %= vertices.size();
    return index;
}
float Intersection::find_default_depth()
{
    int intersect1 = -1, intersect2;
    for (int i = 0; i < vertices.size(); i ++) {
        if (vertices[i].intersect) {
            if (intersect1 == -1) {
                intersect1 = i;
            } else {
                intersect2 = i;
            }
        }
    }

    return find_depth(intersect1, intersect2);
}
float Intersection::find_depth(int start_vertex, int end_vertex)
{
    assert(start_vertex != end_vertex);
    assert_not_equal(vertices[end_vertex].point, vertices[start_vertex].point);
    std::cout << " * find_depth * " << std::endl;
    /*Construct normal */
    /* PERFORMANCE don't know whether we need to normalize */
    glm::vec2 normal = glm::normalize(vertices[end_vertex].point - vertices[start_vertex].point);
    normal = glm::vec2( - normal.y, normal.x );

    /*Construct transformation to align normal to the x-axis*/
    float align_transform_d[4] = {normal.x, -normal.y, normal.y, normal.x};
    glm::mat2 align_transform = make_mat2x2(align_transform_d);

    /* Invariant: the transform makes normal aligned with x-axis, thus start_vertex has maximum y */
    /* So we descend with y-values, to find the max delta x */

    Vertex end(end_vertex, this);
    Vertex backward(start_vertex, this);
    Vertex forward = backward;
    bool current_is_forward = false;

    glm::vec2   forward_p = align_transform * forward->point,
                backward_p = align_transform * backward->point,
                forward_next_p  = align_transform * forward.successive().point,
                backward_next_p = align_transform * backward.preceding().point;

    float biggest_depth = 0;
    do {
        if (forward_next_p.y > backward_next_p.y) {
            ++ forward;
            current_is_forward = true;
            forward_p       = forward_next_p;
            forward_next_p  = align_transform * forward.successive().point;
        } else {
            -- backward;
            current_is_forward = false;
            backward_p      = backward_next_p;
            backward_next_p = align_transform * backward.preceding().point;
        }
        std::cout << "current if forward: " << std::boolalpha << current_is_forward << std::endl;
        std::cout << "\t " << forward.get_index() << ", " << backward.get_index() << std::endl;


        if (current_is_forward) {
            float alpha;
            std::cout << "intersect_horizontal " << backward_p << ", " << backward_next_p << ", " << forward_p.y << std::endl;
            float intersect_x = intersect_horizontal(backward_p, backward_next_p - backward_p, forward_p.y,  alpha);
            std::cout << "alpha=" << alpha << std::endl;
            if (fabs(intersect_x - forward_p.x) > biggest_depth) {
                biggest_depth = fabs(intersect_x - forward_p.x);
            }
        } else {
            float alpha;
            std::cout << "intersect_horizontal " << forward_p << ", " << forward_next_p << ", " << backward_p.y << std::endl;
            float intersect_x = intersect_horizontal(forward_p, forward_next_p - forward_p, backward_p.y,  alpha);
            std::cout << "alpha=" << alpha << std::endl;
            if (fabs(intersect_x - backward_p.x) > biggest_depth) {
                biggest_depth = fabs(intersect_x - backward_p.x);
            }
        
        }
    //float intersect_horizontal(vec2 line_start, vec2 line_direction, float y_constant, float &alpha_out)
    } while (forward != end && backward != end);

    return biggest_depth;
} 
//////////////////////////////////////////////
//              Vertex pointer
//////////////////////////////////////////////

Intersection::Vertex::Vertex(int index, Intersection* parent)
    :index(index % parent->vertices.size()), parent(parent) { }

void Intersection::Vertex::set_index(int val)
{
    index = val;
    if (index < 0) index += parent->vertices.size();
    index = index % parent->vertices.size();
}
HybridVertex& Intersection::Vertex::operator* ()
{
	return parent->vertices[index];
}
HybridVertex* Intersection::Vertex::operator-> ()
{
	return &(parent->vertices[index]);
}
HybridVertex& Intersection::Vertex::successive()
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}
HybridVertex& Intersection::Vertex::preceding()
{
	int i = (index == 0) ? parent->vertices.size() - 1 : index - 1;
	return *(&(parent->vertices[i]));
}
Intersection::Vertex& Intersection::Vertex::operator++ ()
{
	index = (index == parent->vertices.size() - 1) ? 0 : index + 1;
    return *this;
}
Intersection::Vertex& Intersection::Vertex::operator-- ()
{
	index = (index == 0) ? parent->vertices.size() - 1 : index - 1;
    return *this;
}
bool Intersection::Vertex::operator== (Intersection::Vertex& v)
{
    return (index == v.index && parent == v.parent);
}
bool Intersection::Vertex::operator!= (Intersection::Vertex& v)
{
    return (index != v.index || parent != v.parent);
}


///////////
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
/////////// Vertex of Line Strip
// TODO there are a lot of copies of this class for different kinds of geometric things... make template class?

