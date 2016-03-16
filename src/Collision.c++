#include "Polygon.h"
#include "Collision.h"
#include "Geometry.h"
#include "LinkedList.h"
#include "Renderer.h"
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <limits>
#include <utility>

using namespace glm;

LineStrip Intersection::CastInternalShadow(glm::vec2 direction, Polygon* subject)
{
    direction = normalize(direction);
    float align_transform_d[4] = {direction.x, -direction.y, direction.y, direction.x}; //note: column major
    mat2 align_transform = make_mat2x2(align_transform_d);

    /* find the limits of the intersection on the transformed y-axis */
    float y_min = std::numeric_limits<float>::max();
    float y_max = -std::numeric_limits<float>::max();
    Vertex y_min_vertex, y_max_vertex;

    for (int i = 0; i < vertices.size(); i ++) 
    {
        vec2 vertex = align_transform * vertices[i].point;
        if (vertex.y < y_min) {
            y_min = vertex.y;
            y_min_vertex = Vertex(i, this);
        }
        if (vertex.y > y_max) {
            y_max = vertex.y;
            y_max_vertex = Vertex(i, this);
        }
    }
    std::cout << " Y " << y_min << ", " << y_max << std::endl; 

    float alpha;

    EdgePoint start_point;
    start_point.parent = subject;
    EdgePoint end_point;
    end_point.parent = subject;
    bool start_found = false;
    bool end_found = false;
    // (i, j) is the edge indeces, (a, b) is the edge coordinates
    int i = subject->vertices.size() - 1; 
    vec2 a = align_transform * subject->transform(subject->vertices.back());
    for (int j_2 = 0; j_2 < subject->vertices.size() * 2; j_2 ++) {
        int j = j_2 % subject->vertices.size();
        vec2 b = align_transform * subject->transformed(j);
        
        // std::cout << a.y << " < " << y_max << " && " << b.y << " < " << y_max << std::endl;
        if (a.y >= y_max && b.y <= y_max) { // Enter via y_max - start recording the path
            intersect_horizontal(a, b - a, y_max, alpha);
            start_found = true;
            start_point.index = i;
            start_point.alpha = alpha;
        }
        /* std::cout << "Edge (" <<i << ", " << j << ")     - start found: " << std::boolalpha << start_found << std::endl;
        std::cout << a.y << " >= " << y_min << "\t&&\t" << b.y << " <= " << y_min << std::endl; */
        if (a.y >= y_min && b.y <= y_min) { // Leave via y_min
            if (start_found) { // SUCCESS
                intersect_horizontal(a, b - a, y_min, alpha);
                end_found = true;
                end_point.index = i;
                end_point.alpha = alpha;
                break;
            }
        }
        i = j;
        a = b;
    }
    //TODO END NOT FOUND
    /* std::cout << "  start: " << start_point.index << " .. a = " << start_point.alpha << std::endl;
    std::cout << "  end: " << end_point.index << " .. a = " << end_point.alpha << std::endl; */
    assert(start_found);
    assert(end_found);

    /* Create line strip */
	bool CCW = (subject->signedArea() > 0); // Is this intersection conter clockwise?
    LineStrip result(subject);
    EdgePoint p;
    p.parent = subject;

    Polygon::Vertex start(start_point.index, subject), end(end_point.index, subject);

    result.vertices.push_back(start_point);
    if (start_point.index != end_point.index) {

        Polygon::Vertex it = start; 
        ++ it;
        while (true) {
            p.set( it.getIndex(), 0);
            result.vertices.push_back(p);

            if (it == end) break;

            if (CCW)    -- it;
            else        ++ it;
        }
    }
    result.vertices.push_back(end_point);

    return result;
}



Manifold Intersection::manifold(vec2 relative_velocity, Polygon* reference, Polygon* subject, Renderer& renderer)
{
    std::cout << " ---------------- " << std::endl;
    //TODO simplifying assumption: the line strips are both CCW or both CW
    #define v relative_velocity
    float align_transform_d[4] = {v.x, -v.y, v.y, v.x}; //note: column major
    mat2 align_transform = make_mat2x2(align_transform_d);
    mat2 align_inverse = glm::inverse(align_transform);
    std::cout << " -- reference -- " << std::endl;
    LineStrip r_strip = CastInternalShadow(-v, reference);
    std::cout << " -- subject -- " << std::endl;
    LineStrip s_strip = CastInternalShadow(v, subject);
    /*** Task: Find the line of longest distance in direction v, between the two line strips ***/
    /*** So after transformation, find line ... in x direction ***/
    /*** What we do know is that after transformation, the lines start and end at the same y values ***/
    // TODO depending on CCW: We need to know in what direction (+/-) we are moving in the loop
    
    /* Goal: */
    struct ManifoldData {
        float depth = 0;
        // Data for the collision vertex:
        EdgePoint vertex_point; // could also be an EdgePoint but we don't need alpha
        // Data for the collision edge:
        EdgePoint edge_point;

    };
    ManifoldData max_manifold;
    /**/
    LineStrip::Vertex r(r_strip.vertices.size() - 1, &r_strip);
    LineStrip::Vertex s(0, &s_strip);
    vec2 r_vec = align_transform * r->point_t();
    vec2 s_vec = align_transform * s->point_t();
    vec2 r_vec_next, s_vec_next;
#define R false
#define S true
    bool current = R;
    while (true) {
        /* Treat */
        // std::cout << r.getIndex() << ", " << s.getIndex() << std::endl;
        float alpha;
        float intersect_x;
        if (current == R) {
            intersect_x = intersect_horizontal(s_vec, s_vec_next - s_vec, r_vec.y, alpha);
            vec2 projected_point = s_vec + alpha * (s_vec_next - s_vec);
            renderer.addVector(align_inverse * r_vec, align_inverse * (projected_point - r_vec));
            float depth = glm::length(projected_point - r_vec);
            if (depth > max_manifold.depth) {
                max_manifold.depth = depth;
                max_manifold.vertex_point = EdgePoint(r.getIndex(), 0, reference);
                max_manifold.edge_point = EdgePoint(s.getIndex(), alpha, subject);
            }
        } else {
            intersect_x = intersect_horizontal(r_vec, r_vec_next - r_vec, s_vec.y, alpha);
            vec2 projected_point = r_vec + alpha * (r_vec_next - r_vec);
            renderer.addVector(align_inverse * s_vec, align_inverse * (projected_point - s_vec));
            float depth = glm::length(projected_point - s_vec);
            if (depth > max_manifold.depth) {
                max_manifold.depth = depth;
                max_manifold.vertex_point = EdgePoint(s.getIndex(), 0, subject);
                max_manifold.edge_point = EdgePoint(r.getIndex(), alpha, reference);
            }
        }

        /* Advance */
        if (r.atBegin() || s.atEnd()) break;
        r_vec_next = align_transform * (--LineStrip::Vertex(r))->point_t();
        s_vec_next = align_transform * (++LineStrip::Vertex(s))->point_t();

        if (r_vec_next.y > s_vec_next.y) { // We are moving upwards (temporary solution), so pick r because it comes first
            if (!r.atBegin())
                -- r;
            current = R;
        } else {
            if (!s.atEnd())
                ++ s;
            current = S;
        }
        r_vec = align_transform * r->point_t();
        s_vec = align_transform * s->point_t();
    }
    /* std::cout << (max.edge_owner == reference) << std::endl;
    renderer.addDot(align_inverse * max.point->point); */
    Manifold result;
    /*
     * Collision normal and depth are in world coordinates.
     * Collision points are EdgePoints - so local to each polygon
     */
    Polygon::Edge collision_edge(max_manifold.edge_point.index, max_manifold.edge_point.parent);
    vec2 edge_direction = collision_edge.end() - collision_edge.start();
    vec2 collision_normal = vec2(-edge_direction.y, edge_direction.x);

    result.subj_point = max_manifold.vertex_point;
    result.ref_point = max_manifold.edge_point;
    result.normal = collision_normal * max_manifold.depth;

    return result;
}



////////////////////////////
// CALCULATE INTERSECTION //
///////////////////////////

// REMEMBER: q is the small moving one
// TODO MEMORY LEAK? not afaik

std::vector<Intersection> Polygon::ExtractIntersections(Polygon& p, Polygon& q, bool flip_logic)
{
    std::vector<Intersect> intersects = overlaps(p, q); // Let an intersect be an intersection vertex
    bool p_inside_q, q_inside_p;

    /** LINK PARALLEL VERTICES **/
    std::vector<FullIntersect> sorted;
    for (auto it = intersects.begin(); it != intersects.end(); it ++)
    {
        Intersect& i = *it;
        NewVertex *p = new NewVertex(i, OUT);
        NewVertex *q = new NewVertex(i, OUT);
        p->parallel = q; q->parallel = p;
        sorted.push_back(  FullIntersect(&i, p, q)  );
    }

    /** Sort wrt p **/
    //  The FIRST polygon of Intersection is from Polygon p
    std::sort(sorted.begin(), sorted.end(), FullIntersect::lt<Intersect::FIRST>);
    /** Interleave vertices and intersects in p **/
    CircularList<NewVertex> p_vertices;

    Side in_out = p_inside_q = inside(p.transform(p.vertices[0]), q) ^ flip_logic;
    int added_vertex_index = -1;
    int current_index = 0;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge1.getIndex();
        // Add eventual vertices that were skipped
        Vertex vertex(added_vertex_index + 1, &p);
        while (current_index > added_vertex_index) {
            NewVertex *to_add = new NewVertex(vertex, in_out);

            p_vertices.push_back(to_add);

            ++ vertex;
            ++ added_vertex_index;
        }
        // Add the intersection point
        in_out = !in_out;
        it->vert_p->in_out = in_out;
        p_vertices.push_back(it->vert_p);
    }
    // just to test our logic..
    // assert(in_out == inside(p.vertices[0], q));
    // Add rest of vertices..
    if (sorted.size() > 0) {
        for (uint i = sorted.back().i->edge1.getIndex() + 1; i  < p.vertices.size(); i ++)
        {
            NewVertex *to_add = new NewVertex(Vertex(i, &p), in_out);
            p_vertices.push_back(to_add);
        }
    }

    /** Sort wrt q **/
    //  The FIRST polygon of Intersection is from Polygon p
    std::sort(sorted.begin(), sorted.end(), FullIntersect::lt<Intersect::SECOND>);
    /** Interleave vertices and intersects in q **/
    CircularList<NewVertex> q_vertices;

    in_out = q_inside_p = inside(q.transform(q.vertices[0]), p) ^ flip_logic;
    added_vertex_index = -1;
    current_index = 0;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge2.getIndex();
        // Add eventual vertices that were skipped
        Vertex vertex(added_vertex_index + 1, &q);
        while (current_index > added_vertex_index) {
            NewVertex *to_add = new NewVertex(vertex, in_out);
            q_vertices.push_back(to_add);

            ++ vertex;
            ++ added_vertex_index;
        }
        // Add the intersection point
        in_out = !in_out;
        it->vert_q->in_out = in_out;
        q_vertices.push_back(it->vert_q);
    }
    // just to test our logic..
    // assert(in_out == inside(q.vertices[0], p)); // TODO doesn't work when q[0] is inside of p
    // Add rest of vertices..
    if (sorted.size() > 0) {
        for (uint i = sorted.back().i->edge2.getIndex() + 1; i  < q.vertices.size(); i ++)
        {
            NewVertex *to_add = new NewVertex(Vertex(i, &q), in_out);
            q_vertices.push_back(to_add);
        }
    }

    

    std::vector<Intersection> result;
    /** Create Polygons - by traversing q (because we already have the intersection points of q sorted. **/
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        Intersection polygon;
       
        NewVertex* start;
        NewVertex* current = start = it->vert_q;
        // current and start are _intersection_
        
        if (current->processed) continue;

        Direction direction = current->in_out;

        // This is the start of a new polygon.
        do {
            polygon.vertices.push_back(current->vertex);

            // TODO Shouldn't be possible for the new polygon to consist of only two non-adjacent intersection points
            if (current->intersect) {
                assert(current->parallel);
                current->processed = true;
                current->parallel->processed = true;
                current = current->parallel;
                
                direction = current->in_out;
            }

            if (direction == FORTH) {
                assert(current->next);
                current = current->next;
            } else {
                assert(current->prev);
                current = current->prev;
            }
        } while (current != start);

        result.push_back(polygon);
    }
    if (sorted.size() == 0) {
        if (p_inside_q)
            result.push_back( Intersection(p) );
        else if (q_inside_p)
            result.push_back( Intersection(q) );
    }

    return result;
}


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
        std::cout << "Error: edges do not actually overlap (" << edge1.getIndex() << ", " << edge2.getIndex() << ")" << std::endl;
    }
}

template <Intersect::Which which> // TODO separate functions for FIRST and SECOND
bool Intersect::lt(Intersect& i, Intersect& j)
{
    if (which == FIRST) {
        if (i.edge1.getIndex() < j.edge1.getIndex()) {
            return true;
        } else if (i.edge1.getIndex() == j.edge1.getIndex()) {
            return (i.alpha1 < j.alpha1);
        } else {
            return false;
        }
    } else {
        if (i.edge2.getIndex() < j.edge2.getIndex()) {
            return true;
        } else if (i.edge2.getIndex() == j.edge2.getIndex()) {
            return (i.alpha2 < j.alpha2);
        } else {
            return false;
        }
    }
}

template <Intersect::Which which>
bool FullIntersect::lt(FullIntersect& a, FullIntersect& b)
{
    return Intersect::lt<which>(*(a.i), *(b.i));
}


/****************/
/* HybridVertex */
/****************/
HybridVertex::HybridVertex(Intersect& i)
{
    edge1_owner = i.edge1.getParent();
    edge1_index = i.edge1.getIndex();
    edge2_owner = i.edge2.getParent();
    edge2_index = i.edge2.getIndex();
    intersect = true;
    //
    point = i.point;
}
HybridVertex::HybridVertex(Polygon::Vertex v)
{
    owner = v.getParent();
    vertex = v.getIndex();
    intersect = false;
    point = v.transformed();
    // std::cout << "Make HV point " << glm::to_string(point) << std::endl;
}
/****************/
/* Intersection */
/****************/
Intersection::Intersection(Polygon& p)
{
    vertices.resize(p.vertices.size());
    for (int i = 0; i < p.vertices.size(); i ++)
    {
        HybridVertex v;
        v.owner = &p;
        v.vertex = i;
        v.point = p.transform(p.vertices[i]);
        vertices[i] = v;
        // std::cout << p.vertices[i].x << ", " << p.vertices[i].y << std::endl;
    }
    // std::cout << p.vertices.size() << std::endl;
}
void Intersection::appendLinesToVector(std::vector<float> &list)
{
    // Loop through HybridVertices
    std::vector<HybridVertex>::iterator next;
    for (auto it = vertices.begin(); it != vertices.end(); it ++)
    {
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();

        
        glm::vec2 vec_i = it->point;
        glm::vec2 vec_j = next->point;
        list.push_back(vec_i.x);
        list.push_back(vec_i.y);
        list.push_back(vec_j.x);
        list.push_back(vec_j.y);
    }
}

float Intersection::signedArea()
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
	float A = signedArea();
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
//////////////////////////////////////////////
//              Vertex pointer
//////////////////////////////////////////////

Intersection::Vertex::Vertex(int index, Intersection* parent)
    :index(index % parent->vertices.size()), parent(parent) { }

void Intersection::Vertex::setIndex(int val)
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
/////////// Vertex of Line Strip
// TODO there are a lot of copies of this class for different kinds of geometric things... make template class?

LineStrip::Vertex::Vertex(int index, LineStrip* parent)
    :index(index % parent->vertices.size()), parent(parent) { }

void LineStrip::Vertex::setIndex(int val)
{
    index = val;
    if (index < 0) index += parent->vertices.size();
    index = index % parent->vertices.size();
}
glm::vec2 EdgePoint::point()
{
    int next_index = (index + 1) % parent->vertices.size();
    return parent->vertices[index];
}
glm::vec2 EdgePoint::point_t()
{
    int next_index = (index + 1) % parent->vertices.size();
    return parent->transform((1 - alpha) * parent->vertices[index] + alpha * parent->vertices[next_index]);
}
EdgePoint& LineStrip::Vertex::operator* ()
{
	return parent->vertices[index];
}
EdgePoint* LineStrip::Vertex::operator-> ()
{
	return &(parent->vertices[index]);
}
EdgePoint& LineStrip::Vertex::successive()
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}
EdgePoint& LineStrip::Vertex::preceding()
{
	int i = (index == 0) ? parent->vertices.size() - 1 : index - 1;
	return *(&(parent->vertices[i]));
}
LineStrip::Vertex& LineStrip::Vertex::operator++ ()
{
    if (index < parent->vertices.size() - 1) index ++;
    return *this;
}
LineStrip::Vertex& LineStrip::Vertex::operator-- ()
{
    if (index > 0) index --;
    return *this;
}
bool LineStrip::Vertex::operator== (LineStrip::Vertex& v)
{
    return (index == v.index && parent == v.parent);
}
///////////////////
// LineStrip /////
// //////////////
void LineStrip::appendLinesToVector(std::vector<float> &list)
{
    // Loop through HybridVertices
    std::vector<EdgePoint>::iterator prev;
    auto it = vertices.begin();
    prev = it;
    ++ it;
    for (; it != vertices.end(); it ++)
    {
        glm::vec2 vec_i = prev->point_t();
        glm::vec2 vec_j = it->point_t();
        list.push_back(vec_i.x);
        list.push_back(vec_i.y);
        list.push_back(vec_j.x);
        list.push_back(vec_j.y);
        prev = it;
    }
}

///////////////////////////////////////////
// Overlaps: old, maybe not used anymore //
///////////////////////////////////////////
std::vector<Intersect> Polygon::overlaps(Polygon& a, Polygon& b)
{
	// Make an "Influence Area" from a
	// test centroid of b.
	// glm::vec2 centroid_a = a.transform(a.centroid());
	glm::vec2 centroid_b = b.transform(b.centroid());
	float radius_b = b.radius();
	// First, loop through edges (u, v) of a, and create a triangle with the centroid, from which we
	// calculate the barycentric coordinate space
    
    std::vector<Intersect> intersections;

	float distFromEdge;
	float slope, min_val, delta_val;
	glm::vec2 delta;
	bool withinBounds;
	LineSegment edge_a, edge_b;
	for (int i = 0; i < a.numEdges(); i ++)
	{
		edge_a = a.getEdge(i);
        edge_a.first = a.transform(edge_a.first);
        edge_a.second = a.transform(edge_a.second);
		// Get barycentric coordinates with respect to centroid_a
		// glm::vec3 lala = barycentric(centroid_a, edge_a.first, edge_a.second, centroid_b);

		// If absolute value of this is within the radius of b, then there is a possible collision
		
		{	// First, also limit the y or x coordinate based on slope
			delta = edge_a.first - edge_a.second;
			slope = delta.y / delta.x;
			if (fabs(slope) < 1) {	// Limit on x axis
				min_val = std::min(edge_a.first.x, edge_a.second.x) - radius_b;
				delta_val = fabs(delta.x) + radius_b + radius_b;
				withinBounds = centroid_b.x > min_val && centroid_b.x < min_val + delta_val;
			} else {				// Limit on y axis
				min_val = std::min(edge_a.first.y, edge_a.second.y) - radius_b;
				delta_val = fabs(delta.y) + radius_b + radius_b;
				withinBounds = centroid_b.y > min_val && centroid_b.y < min_val + delta_val;
			}
		}
		distFromEdge = distance(centroid_b, edge_a.first, edge_a.second);
		if (fabs(distFromEdge) <= radius_b && withinBounds) {
			// DO DETAILED TEST b vs. this edge (s, t)
			//TODO: eliminate edges where both vertices are outside (using the bary-space of the centroid-triangle)
			//TODO: ... we already know which edges may collide and which not
			for (int j = 0; j < b.numEdges(); j ++)
			{
				edge_b = b.getEdge(j);
                edge_b.first = b.transform(edge_b.first);
                edge_b.second = b.transform(edge_b.second);
				if (intersect(edge_a.first, edge_a.second, edge_b.first, edge_b.second)) {
                    intersections.push_back(Intersect(   Polygon::Edge(i, &a),
                                                            Polygon::Edge(j, &b)));
				}
			}
		}
	}
	return intersections;
}


#if 0
/* IMPORTANT: Take into consideration local coordinate systems */
Intersection Intersection::ExtendInDirection(glm::vec2 v, Polygon* subject)
{
    v = normalize(v);
    float align_transform_d[4] = {v.x, v.y, -v.y, v.x};
    mat2 align_transform = make_mat2x2(align_transform_d);

    /* TODO: transform all polygons. For now, just transform whenever they are used */
    /* find the limits of the intersection on the transformed y-axis */
    float y_min = std::numeric_limits<float>::max();
    float y_max = -std::numeric_limits<float>::max();
    // HybridVertex y_min_vertex;
    // HybridVertex y_max_vertex;
    // int y_min_vertex_index, y_max_vertex_index;
    Vertex y_min_vertex, y_max_vertex;

    for (int i = 0; i < vertices.size(); i ++) 
    {
        vec2 vertex = align_transform * vertices[i].point;
        if (vertex.y < y_min) {
            y_min = vertex.y;
            y_min_vertex = Vertex(i, this);
        }
        if (vertex.y > y_max) {
            y_max = vertex.y;
            y_max_vertex = Vertex(i, this);
        }
    }

    /* Find the portion of the subject which we want to attach to the intersection */
    /* - If clockwise: it has to enter at y-max, and leave at y-min */

    std::vector<HybridVertex> potential_vertices; // list of intersects between y max/min lines and subject

    float alpha;

    HybridVertex start_point;
    HybridVertex end_point;
    bool start_found = false;
    bool end_found = false;
    /* first, if vertex is between a and b, go backwards to find the start point
     * (the first edge we will check next starts at the last vertex - that's why back() )*/
    {
        float first_y = (align_transform * subject->vertices.back()).y;
        if (first_y < y_max && first_y > y_min) {
            int j = 0;
            vec2 b = align_transform * subject->vertices.front();
            for (int i = subject->vertices.size() - 1; i >=0; i --)
            {
                vec2 a = align_transform * subject->vertices[i];
                if (a.y >= y_max && b.y <= y_max) {
                    intersect_horizontal(a, b - a, y_max, alpha);
                    // TODO NEXT: respect alpha in point calculation - need an Intersection::Vertex object...
                    start_found = true;
                    start_point.vertex = i;
                    start_point.owner = subject;
                    start_point.point = subject->getPoint(i, alpha);
                    start_point.alpha = alpha;
                    /* start_point.intersect = false;// TODO can exploit intersect aspect */
                    break;
                }
                j = i;
                b = a;
            }
        }
    }
    // (i, j) is the edge indeces, (a, b) is the edge coordinates
    int i = subject->vertices.size() - 1; 
    vec2 a = align_transform * subject->vertices.back();
    for (int j = 0; j < subject->vertices.size(); j ++) {
        vec2 b = align_transform * subject->vertices[j];
        
        // std::cout << a.y << " < " << y_max << " && " << b.y << " < " << y_max << std::endl;
        if (a.y >= y_max && b.y <= y_max) { // Enter via y_max - start recording the path
            intersect_horizontal(a, b - a, y_max, alpha);
            start_found = true;
            start_point.vertex = i;
            start_point.owner = subject;
            start_point.point = subject->getPoint(i, alpha);
            start_point.alpha = alpha;
        }
        /* std::cout << "Edge (" <<i << ", " << j << ")     - start found: " << std::boolalpha << start_found << std::endl;
        std::cout << a.y << " >= " << y_min << "\t&&\t" << b.y << " <= " << y_min << std::endl; */
        if (a.y >= y_min && b.y <= y_min) { // Leave via y_min
            if (start_found) { // SUCCESS
                intersect_horizontal(a, b - a, y_min, alpha);
                end_found = true;
                end_point.vertex = i;
                end_point.owner = subject;
                end_point.point = subject->getPoint(i, alpha);
                end_point.alpha = alpha;
                break;
            }
        }
        i = j;
        a = b;
    }
    assert(start_found);
    assert(end_found);

	bool CCW = (signedArea() > 0); // Is this intersection conter clockwise?
    /** Make the extended polygon **/
    Intersection result;
    /* First add the left part of the intersection */
// #if 0
    Vertex it = y_min_vertex;
    while (true)
    {
        result.vertices.push_back(*it);
        if (it == y_max_vertex) break;
        if (CCW)    -- it;
        else        ++ it;
    }
// #endif
    /* Then add the rest of the subject polygon */

    // Should not add start vertex if it coincides with the y_max vertex
    // Which (simplified) happens if y_max is one of subject's vertices:
    bool omit_top, omit_bottom;
    if (!y_max_vertex->intersect && y_max_vertex->owner == subject) { // TODO additional test (index is the same) for issue [I]
        std::cout << "omit top - is top" << std::endl;
        omit_top = true;
    }
    // ... or if y_max is an intersect in which one of the owners is subject - this NEEDS an index check
    if (y_max_vertex->intersect) {
        if (y_max_vertex->edge1_owner == subject && y_max_vertex->edge1_index == start_point.vertex) {
            std::cout << "omit top - is intersect" << std::endl;
            omit_top = true;
        }
        if (y_max_vertex->edge2_owner == subject && y_max_vertex->edge2_index == start_point.vertex) {
            std::cout << "omit top - is intersect" << std::endl;
            omit_top = true;
        }
    }
    if (!y_min_vertex->intersect && y_min_vertex->owner == subject) {
        std::cout << "omit bottom - is bottom" << std::endl;
        omit_bottom = true;
    }
    if (y_min_vertex->intersect) {
        if (y_min_vertex->edge1_owner == subject && y_min_vertex->edge1_index == end_point.vertex)
        {
            std::cout << "omit bottom - is intersect" << std::endl;
            omit_bottom = true;
        }
        if (y_min_vertex->edge2_owner == subject && y_min_vertex->edge2_index == end_point.vertex) {
            std::cout << "omit bottom - is intersect" << std::endl;
            omit_bottom = true;
        }
    }
    std::cout << " -------------- " << std::endl;

    /* HybridVertex A;
    A.point = glm::vec2(0, 0);
    result.vertices.push_back(A); */
    /* Add rest of polygon */
    if (!omit_top) {
        result.vertices.push_back(start_point);
    }
    // rest of vertices are normal vertices of the subject polygon
    {
        Polygon::Vertex start(0, subject), end(0, subject);
        if (start_point.edge1_owner == subject) start.setIndex(start_point.edge1_index);
        if (start_point.edge2_owner == subject) start.setIndex(start_point.edge2_index);
        if (end_point.edge1_owner == subject) end.setIndex(end_point.edge1_index);
        if (end_point.edge2_owner == subject) end.setIndex(end_point.edge2_index);
        if (start.getIndex() != end.getIndex()) {
            ++ start; // TODO this needs to change if CCW
            Polygon::Vertex it = start;
            while (true) {
                HybridVertex to_add(it);
                result.vertices.push_back(to_add);
                if (it == end) break;
                if (CCW)    -- it;
                else        ++ it;
            }
        }
    }

    if (!omit_bottom) {
        result.vertices.push_back(end_point);
    }


    //TODO this illustrates a fatal error when the peak of intersection is not the highest point ahead [I]
    /* {
        if (!y_max_vertex.intersect && y_max_vertex.owner == subject) {
            std::cout << "ASSERTION" << std::endl;
            assert(y_max_vertex.vertex == start_point.vertex);
        }
    } */

    return result;
}
#endif
