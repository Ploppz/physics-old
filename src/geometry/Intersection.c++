/** src **/
#include "Polygon.h"
#include "Intersection.h"
#include "geometry.h"
#include "LinkedList.h"
#include "tmp.h"
#include "render/Graphics.h"
#include "debug/debug.h"
#include "error_handling.h"
/****/
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <limits>
#include <utility>

using namespace glm;

extern Graphics *g_graphics;


////////////////////
// Intersection
///////////////////
Intersect::Intersect(Polygon::Edge edge_p, Polygon::Edge edge_q)
{
    owner[P]  = edge_p.get_parent();
    vertex[P] = edge_p.get_index();
    owner[Q]  = edge_p.get_parent();
    vertex[Q] = edge_q.get_index();
    // Find intersection point (last argument is output)
    bool result = intersect(edge_p.start_tr(), edge_p.end_tr(),
            edge_q.start_tr(), edge_q.end_tr(), point, alpha[P], alpha[Q]);
    if (!result) {
        runtime_fatal(Formatter() << "Error: edges do not actually overlap (" << edge_p.get_index() << ", " << edge_q.get_index() << ")");
    }
    /* Fix alphas - shoudn't be too close to 0 (or 1). Because else, intersection_extraction may
     * get the wrong results */
    // TODO... need it here??
}




/****************/
/* HybridVertex */
/****************/
HybridVertex::HybridVertex(Intersect& i)
{
    owner[P] = i.owner[P];
    owner[Q] = i.owner[Q];
    vertex[P] = i.vertex[P];
    vertex[Q] = i.vertex[Q];
    alpha[P] = i.alpha[P];
    alpha[Q] = i.alpha[Q];
    is_active[P] = true;
    is_active[Q] = true;

    point = i.point;
}
HybridVertex::HybridVertex( Polygon::Vertex v, PolygonID pid )
{
    owner[pid] = v.get_parent();
    vertex[pid] = v.get_index();
    alpha[pid] = 0;
    is_active[pid] = true;

    point = v.transformed();
}
EdgePoint HybridVertex::get_edge_point(PolygonID pid)
{
    return EdgePoint(vertex[pid], alpha[pid], owner[pid]);
}


bool HybridVertex::is_intersect()
{
    return is_active[0] && is_active[1];
}

/****************/
/* Intersection */
/****************/
Intersection::Intersection(Polygon& p, PolygonID pid)
{
    vertices.resize(p.num_vertices());
    for (int i = 0; i < p.num_vertices(); i ++)
    {
        HybridVertex v( Polygon::Vertex(i, &p), pid );
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
int Intersection::num_intersects()
{
    int count = 0;
    for (auto it = vertices.begin(); it != vertices.end(); it ++)
    {
        if (it->is_intersect())
            ++ count;
    }
    return count;
}

// TODO: Fundamental flaw with following approach:
// - not always only 2 intersects. may be 4 or any even number.


glm::vec2 Intersection::find_normal_wrt(Polygon* p) {

    int start_vertex = -1, end_vertex;
    for (int i = 0; i < vertices.size(); i ++) {
        if (vertices[i].is_intersect()) {
            if (start_vertex == -1) {
                start_vertex = i;
            } else {
                end_vertex = i;
            }
        }
    }
    assert(vertices[end_vertex].point != vertices[start_vertex].point);
    return find_normal_wrt(p, start_vertex, end_vertex);
}
glm::vec2 Intersection::find_normal_wrt(Polygon* p, int start_vertex, int end_vertex) // or.. find direction of intersect wrt
{
    /* PERFORMANCE don't know whether we need to normalize */
    glm::vec2 normal = glm::normalize(vertices[end_vertex].point - vertices[start_vertex].point);
    normal = glm::vec2( - normal.y, normal.x );

    /*  normal points to the left of (end_vertex - start_vertex)
        found with assumption that Polygon p lies to the left of (end_vertex - start_vertex)
        check if polygon is left of the line .. */
    glm::vec2 next_edge = Vertex(start_vertex, this).successive().point - Vertex(start_vertex, this)->point;
    glm::vec2 prev_edge = Vertex(start_vertex, this).preceding().point - Vertex(start_vertex, this)->point;

    /* Check if next edge is left of the (start_vertex, end_vertex) line */
    if (leftof(next_edge, prev_edge)) {
        if (edge_owner(start_vertex) != p) {
            /* This means that p is at the left op (start_vertex, end_vertex) */
            normal = - normal;
        }
    } else {
        if (edge_owner(start_vertex) == p) {
            /* Same */
            normal = - normal;
        }
    }
    /* Conclusion: normal points _away_ from the boundary of Polygon* p */
    return normal;
}
Polygon* Intersection::edge_owner(int edge_start_index)
{
    HybridVertex v1 = vertices[edge_start_index];
    HybridVertex v2 = vertices[clamp_index(edge_start_index + 1)];

    if ( ! v1.is_intersect() ) {
        return v1.owner[v1.get_active()];
    } else if ( ! v2.is_intersect() ) {
        return v2.owner[v2.get_active()];
    } else {
        /* Find the polygon whose _same single edge_ is present in the two intersections */
       if (v1.vertex[P] == v2.vertex[P] || v1.vertex[P] == v2.vertex[Q]) {
            return v1.owner[P];
        } else {
            /* Fair assumption that now, v1.vertex[Q] must be the same as any of v2's edges */
            assert (v1.vertex[Q] == v2.vertex[P] || v1.vertex[Q] == v2.vertex[Q]);
            return v1.owner[Q];
        }
    }
}
int Intersection::clamp_index(int index)
{
    if (index < 0) index += vertices.size();
    index %= vertices.size();
    return index;
}

struct ManifoldData {
    float depth = 0;
    // Data for the collision vertex:
    EdgePoint vertex_point; // don't need alpha
    // Data for the collision edge:
    EdgePoint edge_point;

};

#if 0
/* CCW sensitive ! */
DepthContact Intersection::get_contact(Polygon* reference, Polygon* subject)
{
    const bool DEBUG = false;
    //assert(num_intersects() == 2); /* Logic: not implemented solution */
    if (DEBUG) std::cout << " NUM INTERSECTS:  " << num_intersects() << std::endl;
    if (num_intersects() > 2) {
        DepthContact dc;
        dc.depth = 0;
        return dc;
    }
    if (DEBUG) std::cout << "GET CONTACT"  << std::endl;
    /* TODO it's critical which way normal points, because we already decided how to iterate the LSSes */
    glm::vec2 normal = find_normal_wrt(reference);
    mat2 align_transform = rotate_coor_system(-normal);

    if (DEBUG) std::cout << "*** Cast shadow on b1 ***" << std::endl;
    LineStripSeries<LEFT> r_strip = cast_shadow_on(reference, normal);
    if (DEBUG) std::cout << "*** Cast shadow on b2 ***" << std::endl;
    LineStripSeries<RIGHT> s_strip = cast_shadow_on(subject, - normal);

    r_strip.set_align_matrix(align_transform);
    s_strip.set_align_matrix(align_transform);

    r_strip.make_y_monotone();
    s_strip.make_y_monotone();
    /*** Task: Find the line of longest distance in direction normal, between the two line strips ***/
    /*** So after transformation, scan lines in x direction ***/
    /*** Note: It's not a guarantee that they will start in the same transformed y-values ***/

#define R false
#define S true
    ManifoldData max_manifold;
    max_manifold.depth = 0;

    bool current;
    bool r_active = false,
         s_active = false;
    LineStripSeries<LEFT>::Vertex<true>     r(&r_strip);
    LineStripSeries<RIGHT>::Vertex<false>   s(&s_strip); 
    LineStripSeries<LEFT>::Vertex<true>     r_next = r;
    LineStripSeries<RIGHT>::Vertex<false>   s_next = s;
    if (r.entry_point().y >= s.entry_point().y) { // r comes before s
        r_active = true;
        ++ r_next;
        current = R;
    }
    if (r.entry_point().y <= s.entry_point().y){
        s_active = true;
        ++ s_next;
        current = S;
    }
    // std::cout << " points y: (" << std::fixed << r.entry_point() << " vs " << s.entry_point() << ")" << std::endl;
    vec2 r_vec, s_vec;
    vec2 r_vec_next, s_vec_next;
    vec2 r_vec_best, s_vec_best;
    bool break_next_loop = false;
    bool max_manifold_found = false; // debug
    while (true) {
        r_vec_next = r_next.entry_point();
        s_vec_next = s_next.entry_point();
        r_vec = r.exit_point();
        s_vec = s.exit_point();
        r_vec_best = r.best_point();
        s_vec_best = s.best_point();
        /* Debug */
        if (DEBUG) {
            std::cout << "Activity: " << r_active << ", " << s_active << std::endl;
            std::cout << "Indices: " << r.get_vertex().get_index() << ", " << s.get_vertex().get_index() << std::endl;
            std::cout << " Y: (" << std::fixed << r.entry_point().y << " vs " << s.entry_point().y << ")" << std::endl;
            std::cout << " Y_next: (" << std::fixed << r_next.entry_point().y << " vs " << s_next.entry_point().y << ")" << std::endl;
        }
        /* Treat */
        if (r_active && s_active) {
            float alpha;
            float intersect_x;
            if (current == R) {
                intersect_x = intersect_horizontal(s_vec, s_vec_next - s_vec, r_vec.y, alpha);
                vec2 projected_point = s_vec + alpha * (s_vec_next - s_vec);


                float depth = glm::length(projected_point - r_vec_best);
                if (DEBUG) std::cout << "DEPTH " << depth << std::endl;
                if (DEBUG) std::cout << "-  Alpha: " << alpha << std::endl;
                if (depth > max_manifold.depth) {
                    if (DEBUG) std::cout << "-  Used " << std::endl;
                    max_manifold.depth = depth;
                    max_manifold.vertex_point = r.to_edge_point(0);
                    max_manifold.edge_point = s.to_edge_point(alpha);
                    max_manifold_found = true;
                }
            } else if (current == S){
                intersect_x = intersect_horizontal(r_vec, r_vec_next - r_vec, s_vec.y, alpha);
                vec2 projected_point = r_vec + alpha * (r_vec_next - r_vec);


                float depth = glm::length(projected_point - s_vec_best);
                if (DEBUG) std::cout << "DEPTH " << depth << std::endl;
                if (DEBUG) std::cout << "-  Alpha: " << alpha << std::endl;
                if (depth > max_manifold.depth) {
                    if (DEBUG) std::cout << "-  Used " << std::endl;
                    max_manifold.depth = depth;
                    max_manifold.vertex_point = s.to_edge_point(0);
                    max_manifold.edge_point = r.to_edge_point(alpha);
                    max_manifold_found = true;
                    if (DEBUG) std::cout << "-  Now: " << max_manifold.vertex_point.alpha << ", " << max_manifold.edge_point.alpha << std::endl;
                }
            }
        }


        //TODO need a new way to break here
        if (break_next_loop) break;
        /** Advance **/


        /* Advance some vertex on current LineStrip - break loop if necessary */
        // std::cout << r_vec_next.y << " vs " << s_vec_next.y << std::endl;
        if (r_vec_next.y > s_vec_next.y) { // We are moving downwards, so pick r because it comes first
            if (DEBUG) std::cout << "Advance R" << std::endl;
            if (r_active)
                ++ r;

            if (!r.at_end()) {
                r_next = r; ++ r_next;
            } else {
                break_next_loop = true;
            }
            current = R;
            r_active = true;
        } else {
            if (DEBUG) std::cout << "Advance S" << std::endl;
            if (s_active) {
                ++ s;
            }
            if (!s.at_end()) {
                s_next = s; ++ s_next;
            } else {
                break_next_loop = true;
            }
            current = S;
            s_active = true;
        }
        // REMEMBER to also process the very end 
    }
    assert(max_manifold_found);
    DepthContact result;
    /*
     * Collision normal is in world coordinates.
     * Collision points are EdgePoints - so local to each polygon
     */
    Polygon::Edge collision_edge(max_manifold.edge_point.index, max_manifold.edge_point.parent);
    vec2 edge_direction = collision_edge.end_tr() - collision_edge.start_tr();
    vec2 collision_normal = normalize(vec2(-edge_direction.y, edge_direction.x));
    if (max_manifold.edge_point.parent->CCW) collision_normal = - collision_normal;

    result.subj_point = max_manifold.vertex_point;
    result.ref_point = max_manifold.edge_point;
    result.normal = collision_normal;
    result.depth = length( result.subj_point.point_t() - result.ref_point.point_t() );

    if (DEBUG) std::cout << "Alphas: " << result.subj_point.alpha << ", " << result.ref_point.alpha << std::endl;
    /* assert( result.subj_point.alpha >= 0 && result.subj_point.alpha <= 1 &&
            result.ref_point.alpha >= 0 && result.ref_point.alpha <= 1); */
    // std::cout << "Index: " << max_manifold.edge_point.index << std::endl;
    // std::cout << "Alpha: " << max_manifold.edge_point.alpha << std::endl;
    vec2 intersection_dir = max_manifold.edge_point.point_t() - max_manifold.vertex_point.point_t();

    g_graphics->extra_line_buffer.add_vector(max_manifold.vertex_point.point_t(), intersection_dir); 

    return result;

}
DepthContact Intersection::get_contact(Polygon* reference, bool ref_outside, Polygon* subject, bool subj_outside)
{
    DebugBegin();
    // normal should be in the direction to move subject out of reference
    vec2 normal = find_normal_wrt(subject);
    const bool DEBUG = false;
    // std::cout << " ---------------- " << std::endl;
    //TODO simplifying assumption: the line strips are both CCW or both CW
#define v normal
#define manifold_debug false
#define R false
#define S true
    float align_transform_d[4] = {v.x, -v.y, v.y, v.x}; //note: column major
    mat2 align_transform = make_mat2x2(align_transform_d);
    mat2 align_inverse = glm::inverse(align_transform);
    /* Extract linstrips */
    float factor;
    if (ref_outside) factor = -1;
    else factor = 1;
    LineStripSeries<LEFT> r_strip = cast_internal_shadow(normal * factor, reference);
    if (subj_outside) factor = -1;
    else factor = 1;
    LineStripSeries<RIGHT> s_strip = cast_internal_shadow(-normal * factor, subject);
    r_strip.set_align_matrix(align_transform);
    s_strip.set_align_matrix(align_transform);

    r_strip.make_y_monotone();
    s_strip.make_y_monotone();
    /**/
    /*** Task: Find the line of longest distance in direction v, between the two line strips ***/
    /*** So after transformation, scan lines in x direction ***/
    /*** (After transformation, the lines start and end at the same y values) ***/
    // TODO depending on CCW: We need to know in what direction (+/-) we are moving in the loop
    
    /* Goal: */
    ManifoldData max_manifold;
    max_manifold.depth = 0;
    /**/

    /* We have a series of LineStrips to iterate through.. */

    if (DEBUG) {
        /* std::cout << "Reference: " << r_it->start.index << std::endl;
        std::cout << "Subject: " << s_it->start.index << std::endl;  */
    }
    LineStripSeries<LEFT>::Vertex<true> r(&r_strip); // TODO set to end
    LineStripSeries<RIGHT>::Vertex<false> s(&s_strip); // TODO set to start
    if (DEBUG) dout << "CHECK 1 " <<s.parent->get_parent() << newl;
    if (DEBUG) dout << "CHECK 2 " <<r.parent->get_parent() << newl;
    LineStripSeries<LEFT>::Vertex<true> r_next = r; ++ r_next;
    LineStripSeries<RIGHT>::Vertex<false> s_next = s; ++ s_next;
    vec2 r_vec, s_vec;
    vec2 r_vec_next, s_vec_next;
    vec2 r_vec_best, s_vec_best;
    bool current = R;
    bool break_next_loop = false;
    bool max_manifold_found = false; // debug
    while (true) {
        r_vec_next = r_next.entry_point();
        s_vec_next = s_next.entry_point();
        r_vec = r.exit_point();
        s_vec = s.exit_point();
        r_vec_best = r.best_point();
        s_vec_best = s.best_point();
        /* Debug */
        dout << ".. get_contact iter " << newl;
        /* Treat */
        float alpha;
        float intersect_x;
        if (current == R) {
            intersect_x = intersect_horizontal(s_vec, s_vec_next - s_vec, r_vec.y, alpha);
            vec2 projected_point = s_vec + alpha * (s_vec_next - s_vec);

            if (manifold_debug) g_graphics->extra_line_buffer.add_vector(align_inverse * r_vec, align_inverse * (projected_point - r_vec));

            float depth = glm::length(projected_point - r_vec_best);
            dout << "DEPTH " << depth << newl;
            if (depth > max_manifold.depth) {
                max_manifold.depth = depth;
                max_manifold.vertex_point = r.to_edge_point(0);
                max_manifold.edge_point = s.to_edge_point(alpha);
                max_manifold_found = true;
            }
        } else {
            intersect_x = intersect_horizontal(r_vec, r_vec_next - r_vec, s_vec.y, alpha);
            vec2 projected_point = r_vec + alpha * (r_vec_next - r_vec);

            if (manifold_debug) g_graphics->extra_line_buffer.add_vector(align_inverse * s_vec, align_inverse * (projected_point - s_vec));

            float depth = glm::length(projected_point - s_vec_best);
            dout << "DEPTH " << depth << newl;
            if (depth > max_manifold.depth) {
                max_manifold.depth = depth;
                max_manifold.vertex_point = s.to_edge_point(0);
                max_manifold.edge_point = r.to_edge_point(alpha);
                max_manifold_found = true;
            }
        }


        //TODO need a new way to break here
        if (break_next_loop) break;
        /** Advance **/


        /* Advance some vertex on current LineStrip - break loop if necessary */
        if (r_vec_next.y > s_vec_next.y) { // We are moving downwards (temporary solution), so pick r because it comes first
            ++ r;
            if (!r.at_end()) { // TODO do we need these if's ?
                r_next = r; ++ r_next;
            } else {
                break_next_loop = true;
            }
            current = R;
        } else {
            ++ s;
            if (!s.at_end()) {
                s_next = s; ++ s_next;
            } else {
                break_next_loop = true;
            }
            current = S;
        }
        // REMEMBER to also process the very end 
    }
    assert(max_manifold_found);
    DepthContact result;
    /*
     * Collision normal is in world coordinates.
     * Collision points are EdgePoints - so local to each polygon
     */
    dout << "Parent comparision: " << subject << ", " << reference << ", " <<max_manifold.edge_point.parent << newl;
    Polygon::Edge collision_edge(max_manifold.edge_point.index, max_manifold.edge_point.parent);
    vec2 edge_direction = collision_edge.end_tr() - collision_edge.start_tr();
    vec2 collision_normal = normalize(vec2(-edge_direction.y, edge_direction.x));

    result.subj_point = max_manifold.vertex_point;
    result.ref_point = max_manifold.edge_point;
    result.normal = collision_normal;
    result.depth = length( result.subj_point.point_t() - result.ref_point.point_t() );

    // std::cout << "Index: " << max_manifold.edge_point.index << std::endl;
    // std::cout << "Alpha: " << max_manifold.edge_point.alpha << std::endl;
    vec2 intersection_dir = max_manifold.edge_point.point_t() - max_manifold.vertex_point.point_t();

    g_graphics->extra_line_buffer.add_vector(max_manifold.vertex_point.point_t(), intersection_dir); 

    return result;
#undef v
} 
#endif

bool index_is_next_given_range(int index, int next_index, int range)
{
    int anticipated_next_index = (index + 1) % range;
    return anticipated_next_index == next_index;
}
PolygonID Intersection::find_parent_of_edge(Vertex v1, Vertex v2) // TODO: Only need one Vertex.
{
    if ( !v1->is_intersect() )
        return v1->get_active();
    else if ( !v2->is_intersect() )
        return v2->get_active();
    else
        if ( index_is_next_given_range(v1->vertex[P], v2->vertex[P], v1->owner[P]->num_vertices()) ) {
            return P;
        } else {
            return Q;
        }
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



std::ostream& operator<< (std::ostream& out, Intersection& I)
{
    out << "Intersection: " << std::endl;
    for (int i = 0; i < I.vertices.size(); i ++)
    {
        HybridVertex v = I.vertices[i];
        out << "\t" << i << ": ";
        out << "(" << v.vertex[P] << ", " << v.alpha[1] << ")";
        if (v.is_intersect()) {
            out << " intersects (" << v.vertex[Q] << ", " << v.alpha[2] << ")";
        }
        out << std::endl;
    }
    return out;
}
void Intersection::append_lines_to_vector(std::vector<float>& buffer)
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
        buffer.push_back(vec_i.x);
        buffer.push_back(vec_i.y);
        buffer.push_back(vec_j.x);
        buffer.push_back(vec_j.y);
    }
}
glm::vec2 calculate_point(HybridVertex v) {
    Polygon::Edge edge(v.vertex[P], v.owner[P]);
    return (1.f - v.alpha[P]) * edge.start_tr() + v.alpha[P] * edge.end_tr();
}
void Intersection::append_lines_to_vector2(std::vector<float>& buffer)
{
    // Loop through HybridVertices
    std::vector<HybridVertex>::iterator next;
    for (auto it = vertices.begin(); it != vertices.end(); it ++)
    {
		next = it; next ++;
		if (next == vertices.end())
            next = vertices.begin();

        
        std::cout << "Start: " << it->vertex << ", " << it->alpha << std::endl;
        glm::vec2 vec_i = calculate_point(*it);
        glm::vec2 vec_j = calculate_point(*next);
        // vec_i = it->point;
        // vec_j = next->point;
        buffer.push_back(vec_i.x);
        buffer.push_back(vec_i.y);
        buffer.push_back(vec_j.x);
        buffer.push_back(vec_j.y);
    }
}
