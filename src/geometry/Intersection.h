#pragma once
#include "Polygon.h"
#include "constants.h"
#include "EdgePoint.h"
#include "algorithm/Contact.h"

class LineStrip;

struct Intersect;

struct HybridVertex;
struct Intersection;
enum PolygonID { P=0, Q=1 };

// External:
class Renderer;

struct Intersect
{
    Intersect(Polygon::Edge edge_p, Polygon::Edge edge_q);
    Polygon* owner[2];
    int vertex[2];
    float alpha[2];
#if 0 // old
    Polygon::Edge edge1, edge2;
    float alpha1, alpha2;
#endif
    glm::vec2 point;

    template <PolygonID which>
    static bool lt(Intersect& i, Intersect& j);
};

//////////////////
// INTERSECTION //
//////////////////

// HybridVertex can be an intersect, or can have info about only one edge point
// It may contain information about polygon `P` and/or polygon `Q`.

struct HybridVertex
{
    HybridVertex(){} 
    HybridVertex( Intersect& i );
    HybridVertex( Polygon::Vertex v, PolygonID pid );
    // queries for ease //
    EdgePoint get_edge_point(PolygonID pid);
    bool is_intersect();
    PolygonID get_active(); // TODO implement
    PolygonID get_polygon_id(Polygon* ptr);
    // data //
    Polygon* owner[2];
    int vertex[2];
    float alpha[2];
    bool in_out[2];

    bool is_active[2];

    glm::vec2 point; // only a sort of cache - can be calculated any time
#if 0
    //// old:
    union {
        Polygon* owner,
               * edge1_owner;
    };
    union {
        int vertex,
            edge1_index;
    };
    union {
        float alpha,
              alpha1;
    };

    // For intersects only:
    bool intersect:1;
    bool edge1_in_out:1;
    bool edge2_in_out:1;
    Polygon* edge2_owner;
    int edge2_index;
    float alpha2;
#endif
};


struct Intersection
{
    class Vertex;

    std::vector<HybridVertex> vertices;

    Intersection() {};
    Intersection(Polygon& p, PolygonID pid); // 'copy' p

    // Query //
    PolygonID find_parent_of_edge(Vertex v1, Vertex v2);

    float signed_area();
    glm::vec2 centroid();
    glm::vec2 get_point(int vertex_number, float alpha);

    /* Direction where the polygon is on the Intersection */
    glm::vec2 find_normal_wrt(Polygon* polygon, int start_vertex, int end_vertex);
    glm::vec2 find_normal_wrt(Polygon* polygon);

    /* get_contact:
     * - construct normal from two intersects
     * - project out of each other in that direction
     * - "subj" = polygon of vertex, "ref" = polygon of edge in collision
     * - normal points out of "ref"
     * - returns depth=0 if erroneous/unimplemented case
     */
    /* new */
    DepthContact get_contact(Polygon* reference, Polygon* subject);
    /* old */
    DepthContact get_contact(Polygon* reference, bool ref_outside, Polygon* subject, bool subj_outside);

    Polygon* edge_owner(int edge_start_index);

    LineStrip cast_shadow_on(Polygon* polygon, glm::vec2 direction);
    /* Old: */
    LineStrip cast_internal_shadow(glm::vec2 direction, Polygon* subject);

    void append_lines_to_vector(std::vector<float>& buffer);
    void append_lines_to_vector2(std::vector<float>& buffer);
    
    int num_intersects();

    class Vertex {
    public:
        Vertex(): index(0), parent(0) {};
        Vertex(int index, Intersection *parent);
        int get_index() { return index; }
        void set_index(int val);
        Intersection* get_parent() { return parent; }
        
        HybridVertex & operator* ();
        HybridVertex * operator-> ();

        // Changes what vertex is pointed to
        Vertex& operator++ ();
        Vertex& operator-- ();
        //
        bool operator== (Vertex& v);
        bool operator!= (Vertex& v);

        HybridVertex & preceding();
        HybridVertex & successive();

    private:
        int index;
        Intersection *parent;
    };

    private:
    int clamp_index(int index);
    float find_depth(int start_vertex, int end_vertex);

    EdgePoint interpolate(Vertex v1, Vertex v2, float alpha);

};
std::ostream& operator<< (std::ostream&, Intersection&);




/*** Templates ***/

template <PolygonID which> // TODO separate functions for FIRST and SECOND
bool Intersect::lt(Intersect& i, Intersect& j)
{
    if (i.vertex[which] < j.vertex[which]) {
        return true;
    } else if (i.vertex[which] == j.vertex[which]) {
        return (i.alpha[which] < j.alpha[which]);
    } else {
        return false;
    }
}
