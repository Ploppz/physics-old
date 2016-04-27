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

glm::vec2 Intersection::find_normal_wrt(Polygon* p, int start_vertex, int end_vertex) // or.. find direction of intersect wrt
{
    /* PERFORMANCE don't know whether we need to normalize */
    glm::vec2 normal = glm::normalize(vertices[end_vertex].point - vertices[start_vertex].point);
    normal = glm::vec2( - normal.y, normal.x );

    /*  normal points to the left 
        found with assumption that Polygon p lies to the left of (end_vertex - start_vertex)
        check if polygon is left of the line .. */
    glm::vec2 next_edge = Vertex(start_vertex, this).successive().point - Vertex(start_vertex, this)->point;
    glm::vec2 prev_edge = Vertex(start_vertex, this).preceding().point - Vertex(start_vertex, this)->point;
    if (leftof(next_edge, prev_edge)) {
        if (edge_owner(start_vertex) == p) {
            normal = - normal;
        }
    } else {
        if (edge_owner(start_vertex) != p) {
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
IntersectionContact Intersection::get_contact(Polygon* subject)
{
    IntersectionContact contact;
    int start_vertex = -1, end_vertex;
    for (int i = 0; i < vertices.size(); i ++) {
        if (vertices[i].intersect) {
            if (start_vertex == -1) {
                start_vertex = i;
            } else {
                end_vertex = i;
            }
        }
    }
    if (vertices[end_vertex].point == vertices[start_vertex].point) {
        exit(0);
    }

    /* Construct normal */
    /* PERFORMANCE don't know whether we need to normalize */
    contact.normal = find_normal_wrt(subject, start_vertex, end_vertex);

    /* PERFORMANCE this is used in calculation of normal... also: don't know whether it should be normalized. */
    glm::vec2 tangent = glm::normalize(vertices[end_vertex].point - vertices[start_vertex].point);
    /*Construct transformation to align normal to the x-axis*/
    /* This transform makes start_vertex leftmost and end_vertex rightmost */
    float align_transform_d[4] = {tangent.x, -tangent.y, tangent.y, tangent.x};
    glm::mat2 align_transform = make_mat2x2(align_transform_d);

    /* Invariant: the transform makes normal aligned with x-axis, thus start_vertex has maximum y */
    /* So we descend with y-values, to find the max delta x */

    Vertex end(end_vertex, this);
    Vertex backward(start_vertex, this);
    Vertex forward = backward;
    bool current_is_forward = false;
    bool forward_is_subject = ( edge_owner(forward.get_index()) == subject );

    glm::vec2   forward_p = align_transform * forward->point,
                backward_p = align_transform * backward->point,
                forward_next_p  = align_transform * forward.successive().point,
                backward_next_p = align_transform * backward.preceding().point;

    std::cout << "Forward x: " << forward_p.x << ", " << forward_next_p.x << std::endl;
    std::cout << "Backward x: " << backward_p.x << ", " << backward_next_p.x << std::endl;
    contact.depth = 0;
    do {
        if (forward_next_p.x < backward_next_p.x) {
            std::cout << "Advance forward " << std::endl;
            ++ forward;
            current_is_forward = true;
            forward_p       = forward_next_p;
            forward_next_p  = align_transform * forward.successive().point;
        } else {
            std::cout << "Advance backward " << std::endl;
            -- backward;
            current_is_forward = false;
            backward_p      = backward_next_p;
            backward_next_p = align_transform * backward.preceding().point;
        }
        std::cout << "Forward x: " << forward_p.x << ", " << forward_next_p.x << std::endl;
        std::cout << "Backward x: " << backward_p.x << ", " << backward_next_p.x << std::endl;

        if (current_is_forward) {
            float alpha;
            float intersect_y = intersect_vertical(backward_p, backward_next_p - backward_p, forward_p.x,  alpha);
            std::cout << ":alpha = " << alpha << std::endl;
            // assert(alpha >= 0 && alpha <= 1);
            if (alpha < 0 || alpha > 0) {
                continue;
            }
            if (fabs(intersect_y - forward_p.y) > contact.depth) {
                contact.depth = fabs(intersect_y - forward_p.y);
                //
                contact.subj_point = EdgePoint(forward->vertex, forward->alpha, forward->owner);
                Vertex backward_minus_one = backward; -- backward_minus_one;
                contact.ref_point = interpolate(backward_minus_one, backward, alpha);
                //
                if (forward->owner != subject)
                    std::swap(contact.subj_point, contact.ref_point);
            }
        } else {
            float alpha;
            float intersect_y = intersect_vertical(forward_p, forward_next_p - forward_p, backward_p.x,  alpha);
            std::cout << "alpha = " << alpha << std::endl;
            // assert(alpha >= 0 && alpha <= 1);
            if (alpha < 0 || alpha > 0) {
                continue;
            }
            if (fabs(intersect_y - backward_p.y) > contact.depth) {
                contact.depth = fabs(intersect_y - backward_p.y);
                //
                contact.subj_point = EdgePoint(backward->vertex, backward->alpha, backward->owner);
                Vertex forward_minus_one = forward; -- forward_minus_one;
                contact.ref_point = interpolate(forward_minus_one, forward, alpha);
                //
                if (backward->owner != subject)
                    std::swap(contact.subj_point, contact.ref_point);
            }
        
        }
    //float intersect_horizontal(vec2 line_start, vec2 line_direction, float y_constant, float &alpha_out)
    } while (forward != end && backward != end);

    std::cout << "ref edgepoint: " <<  contact.ref_point.alpha << std::endl;
    std::cout << "subj edgepoint: " << contact.subj_point.alpha << std::endl;
    return contact;
} 

bool index_is_next_given_range(int index, int next_index, int range)
{
    int anticipated_next_index = (index + 1) % range;
    return anticipated_next_index == next_index;
}
Polygon* Intersection::find_parent_of_intersection_edge(Vertex v1, Vertex v2)
{
    /* Assumption: Both vertices are intersects */
    int v1_first_owner_index1 = v1->edge1_index;
    int v1_first_owner_index2;
    if (v1->edge1_owner == v2->edge1_owner) {
        if (index_is_next_given_range(v1->edge1_index, v2->edge1_index, v1->edge1_owner->vertices.size()) )
            return v1->edge1_owner;
        else // condition must be true for the other polygon
            return v1->edge2_owner;
    } else {
        if (index_is_next_given_range(v1->edge1_index, v2->edge2_index, v1->edge1_owner->vertices.size()) )
            return v1->edge1_owner;
        else // condition must be true for the other polygon
            return v1->edge2_owner;
    }
}
EdgePoint Intersection::interpolate(Vertex v1, Vertex v2, float alpha)
{
    EdgePoint e1, e2;
    if (v1->intersect) {
        if (v2->intersect) {
            Polygon* edge_parent = find_parent_of_intersection_edge(v1, v2);
            if (v1->edge1_owner == edge_parent) {
                e1 = EdgePoint(v1->edge1_index, v1->alpha1, edge_parent);
            } else {
                e1 = EdgePoint(v1->edge2_index, v1->alpha2, edge_parent);
            }
            if (v2->edge1_owner == edge_parent) {
                e2 = EdgePoint(v2->edge1_index, v2->alpha1, edge_parent);
            } else {
                e2 = EdgePoint(v2->edge2_index, v2->alpha2, edge_parent);
            }
        } else {
            e2 = EdgePoint(v2->vertex, v2->alpha, v2->owner);

            if (v1->edge1_owner == v2->owner)
                e1 = EdgePoint(v1->edge1_index, v1->alpha1, v2->owner);
            else // v2->edge2_owner == v2->owner
                e1 = EdgePoint(v1->edge2_index, v1->alpha2, v2->owner);
        }
    } else {
        e1 = EdgePoint(v1->vertex, v1->alpha, v1->owner);
        if (v2->intersect) {
            if (v2->edge1_owner == v1->owner)
                e2 = EdgePoint(v2->edge1_index, v2->alpha1, v1->owner);
            else // v2->edge2_owner == v1->owner
                e2 = EdgePoint(v2->edge2_index, v2->alpha2, v1->owner);
        } else {
            e2 = EdgePoint(v2->vertex, v2->alpha, v2->owner);
        }
    }
    assert(e1.parent == e2.parent);
    assert(e1.index == e2.index || e2.alpha == 0);

    /* Now, actually interpolate the EdgePoints ... */
    float alpha1 = e1.alpha;
    float alpha2 = e2.alpha;
    if (e2.index == e1.index) alpha2 = 1;
    float new_alpha = alpha1 + alpha * (alpha2 - alpha1);
    return EdgePoint(e1.index, new_alpha, e1.parent);
    
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

