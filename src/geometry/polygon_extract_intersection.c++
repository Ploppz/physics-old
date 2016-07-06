#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <cstdio> // DELETE 
#include <iomanip> // DELETE
/* src */
#include "Polygon.h"
#include "geometry.h"
#include "constants.h"
#include "Intersection.h"
#include "LinkedList.h"
#include "debug/debug.h"


#include "render/Renderer.h"
#include "debug/StatisticsCollection.h"
extern Renderer* g_renderer;
extern StatisticsCollection* g_statistics;

/** Structs used for algorithm */
struct NewVertex {
    NewVertex(Intersect& i, Side in_out) {
        point = i.point;
        intersect = true;
        processed = false;
        this->in_out = in_out;

        vertex = HybridVertex(i);
        // std::cout << "Creating NewVertex 1, alphas: " << vertex.alpha << ", " << vertex.alpha2 << ";" << std::endl;
    }
    NewVertex(Polygon::Vertex vert, Side in_out) {
        point = vert.transformed();
        intersect = false;
        this->in_out = in_out;

        vertex = HybridVertex(vert);
        // std::cout << "Creating NewVertex 2, alphas: " << vertex.alpha << ", " << vertex.alpha2 << ";" << std::endl;
    }

    glm::vec2 point;
    Side in_out; // Entry or exit to other polygon?
    /* Only for intersects: */
    bool intersect; // Intersection point or normal vertex
    bool processed;

    /* Needed for  linked list */
    NewVertex *prev, *next, *parallel; // Linked list, and link to evt. same intersection point in other polygon

    /* Needed for useful intersection extraction */
    HybridVertex vertex; // it could be possible for the algorithm itself to keep track of this data
                        // e.g. polygon and edge numbers
};

struct FullIntersect { // Needed to sort the NewVertices
    FullIntersect(Intersect *i, NewVertex *p, NewVertex *q): i(i), vert_p(p), vert_q(q){}
    Intersect* i;
    NewVertex* vert_p;
    NewVertex* vert_q;
    template <Intersect::Which which> // sort wrt p or q?
    static bool lt(FullIntersect& a, FullIntersect& b);
};
////////////////////////////
// CALCULATE INTERSECTION //
///////////////////////////

std::vector<Intersection> Polygon::extract_intersections(Polygon& p, Polygon& q, bool p_inside_out, bool q_inside_out)
{
    DebugBeginC(false);
    const bool DEBUG = false;
    std::vector<Intersect> intersects = find_intersects(p, q); // Let an intersect be an intersection vertex
    dout << "Number of intersects: " << intersects.size() << newl;
    if (intersects.size() % 2 != 0 || intersects.size() == 0) {
        return std::vector<Intersection> {};
    }

    if (DEBUG)
    { // test
        std::cout << "INTERSECTS IN : " << std::endl;
        for (auto it = intersects.begin(); it != intersects.end(); it ++) 
        {
            std::cout << it->point << std::endl;
            std::cout << "\t " << it->alpha1 << ", " << it->alpha2 << std::endl;
        }
        std::cout << "END OF LIST " << std::endl;
    }
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


    Side in_out = p_inside_q = inside_stable(p.transformed(0), q) ^ q_inside_out;
    int added_vertex_index = -1;
    int current_index = 0;
    
    dout << "(p)Initial in-out " << std::boolalpha << in_out << newl;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge1.get_index();
        dout << "-- iteration, with vertex index " << current_index << newl;
        // Add eventual vertices that were skipped
        Vertex vertex(added_vertex_index + 1, &p);
        while (current_index > added_vertex_index) {
            NewVertex *to_add = new NewVertex(vertex, in_out);
            dout << " - add normal vertex " << vertex.get_index() << ", in_out = " << in_out << newl;

            p_vertices.push_back(to_add);

            ++ vertex;
            ++ added_vertex_index;
        }
        // Add the intersection point
        in_out = !in_out;
        it->vert_p->in_out = in_out;
        p_vertices.push_back(it->vert_p);
        dout << " - add intersect... (?)" << newl;

        if (in_out)
            dout << " - in" << newl;
        else
            dout << " - out" << newl;
    }
    // just to test our logic..
    assert(in_out == (inside_stable(p.transformed(0), q) ^ q_inside_out));
    // Add rest of vertices..
    if (sorted.size() > 0) {
        for (uint i = sorted.back().i->edge1.get_index() + 1; i  < p.vertices.size(); i ++)
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

    in_out = q_inside_p = inside_stable(q.transformed(0), p) ^ p_inside_out;
    added_vertex_index = -1;
    current_index = 0;

    dout << "(q)Initial in-out " << std::boolalpha << in_out << newl;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge2.get_index();
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
        dout << " - " << in_out << newl;
    }
    // just to test our logic..
    // assert(in_out == inside(q.vertices[0], p)); // TODO doesn't work when q[0] is inside of p
    assert(in_out == (inside_stable(q.transformed(0), p) ^ p_inside_out));
    // Add rest of vertices..
    if (sorted.size() > 0) {
        for (uint i = sorted.back().i->edge2.get_index() + 1; i  < q.vertices.size(); i ++)
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

        Direction direction = (Direction)current->in_out;

        // This is the start of a new polygon.
        bool invalid = false;
        glm::vec2 prev_added_point(0);
        do {
            if (length_squared(current->vertex.point - prev_added_point) < 0.0001f) {
                invalid = true;
            }
            polygon.vertices.push_back(current->vertex);
            prev_added_point = current->vertex.point;

            // (?) TODO Shouldn't be possible for the new polygon to consist of only two non-adjacent intersection points
            if (current->intersect) {
                assert(current->parallel);
                current->processed = true;
                current->parallel->processed = true;
                current = current->parallel;
                
                direction = (Direction)current->in_out;
            }

            if (direction == FORTH) {
                assert(current->next);
                current = current->next;
            } else {
                assert(current->prev);
                current = current->prev;
            }
        } while (current != start);

        if (! invalid) {
            result.push_back(polygon);
        }
    }
    if (sorted.size() == 0) {
        if (p_inside_q) {
            result.push_back( Intersection(p) );
            dout << "p insize q" << newl;
        }
        else if (q_inside_p) {
            result.push_back( Intersection(q) );
            dout << "q insize p" << newl;
        }
    }
    if (true)
    {
        dout << "Found this: " << newl;
        for (auto it = result.begin(); it != result.end(); it ++)
        {
            dout << "# vertices: " << it->vertices.size() << newl;

        }
    }

    return result;
}


std::vector<Intersect> Polygon::find_intersects(Polygon& a, Polygon& b)
{
    std::vector<Intersect> intersections;
#if 0
    LineSegment edge_a, edge_b;
    ////
    edge_a.end = a.transformed(a.num_edges() - 1);
    for (int i = 0; i < a.num_edges(); i ++) {
        edge_a.start = edge_a.end;
        edge_a.end = a.transformed(i);
        ////
        edge_b.end = b.transformed(b.num_edges() - 1);
        for (int j = 0; j < b.num_edges(); j ++) {
            edge_b.start = edge_b.end;
            edge_b.end = b.transformed(j);
            // std::cout << "B: " << edge_b.start << " ; " << edge_b.end << std::endl;
            ////
            glm::vec2 point_of_intersection; float alpha1; float alpha2;
            if (intersect(edge_a.start, edge_a.end, edge_b.start, edge_b.end, point_of_intersection, alpha1, alpha2)) {
                intersections.push_back(Intersect(   Polygon::Edge(i-1, &a),
                                                     Polygon::Edge(j-1, &b)));
            }
        }
    }
    return intersections;

#else

    DebugBegin();
	// Make an "Influence Area" from a
	// test centroid of b.
	// glm::vec2 centroid_a = a.transform(a.centroid());
	glm::vec2 centroid_b = b.transform(b.centroid());
	float radius_b = b.radius();
	// First, loop through edges (u, v) of a, and create a triangle with the centroid, from which we
	// calculate the barycentric coordinate space
    

	float dist_from_edge;
	float min_val, delta_val;
	glm::vec2 delta;
	bool within_bounds;
	LineSegment edge_a, edge_b;
	for (int i = 0; i < a.num_edges(); i ++)
	{
		edge_a = a.get_edge(i);
        edge_a.start = a.transform(edge_a.start);
        edge_a.end = a.transform(edge_a.end);
		// Get barycentric coordinates with respect to centroid_a
		// glm::vec3 lala = barycentric(centroid_a, edge_a.start, edge_a.end, centroid_b);

		// If absolute value of this is within the radius of b, then there is a possible collision
		
		{	// First, also limit the y or x coordinate based on slope
			delta = edge_a.start - edge_a.end;
			if (fabs(delta.x) > fabs(delta.y)) {	// Limit on x axis
				min_val = std::min(edge_a.start.x, edge_a.end.x) - radius_b;
				delta_val = fabs(delta.x) + radius_b + radius_b;
				within_bounds = centroid_b.x > min_val && centroid_b.x < min_val + delta_val;
			} else {				// Limit on y axis
				min_val = std::min(edge_a.start.y, edge_a.end.y) - radius_b;
				delta_val = fabs(delta.y) + radius_b + radius_b;
				within_bounds = centroid_b.y > min_val && centroid_b.y < min_val + delta_val;
			}
		}
		dist_from_edge = distance_line_segment(centroid_b, edge_a.start, edge_a.end);
		if (fabs(dist_from_edge) <= radius_b && within_bounds) {
			// DO DETAILED TEST b vs. this edge (s, t)
			//TODO: eliminate edges where both vertices are outside (using the bary-space of the centroid-triangle)
			//TODO: ... we already know which edges may collide and which not
			for (int j = 0; j < b.num_edges(); j ++)
			{
				edge_b = b.get_edge(j);
                edge_b.start = b.transform(edge_b.start);
                edge_b.end = b.transform(edge_b.end);

                glm::vec2 point_of_intersection;
                float alpha1;
                float alpha2;
				if (intersect(edge_a.start, edge_a.end, edge_b.start, edge_b.end, point_of_intersection, alpha1, alpha2)) {
                    intersections.push_back(Intersect(   Polygon::Edge(i, &a),
                                                            Polygon::Edge(j, &b)));
				}
			}
		}
	}
	return intersections;
#endif
}




template <Intersect::Which which>
bool FullIntersect::lt(FullIntersect& a, FullIntersect& b)
{
    return Intersect::lt<which>(*(a.i), *(b.i));
}
