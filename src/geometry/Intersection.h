#pragma once
#include "Polygon.h"
#include "constants.h"
#include "EdgePoint.h"

class LineStrip;

struct Intersect;

struct HybridVertex;
struct Intersection;

// External:
class Renderer;

struct Intersect
{
    Intersect(Polygon::Edge edge1, Polygon::Edge edge2);
    Polygon::Edge edge1, edge2;
    glm::vec2 point;
    float alpha1, alpha2; // How far along the edges the point resides

    // Lessthan function -- choose which to compare, edge1 or edge2
    enum Which { FIRST, SECOND};
    template <Which which>
    static bool lt(Intersect& i, Intersect& j);
};

//////////////////
// INTERSECTION //
//////////////////

/* HybridVertex takes two forms:
    (owner, vertex, alpha) for regular vertices, and (edge1_owner, edge1_index, alpha1, edge2_owner, edge2_index, alpha2) for intersects
    (coor is for all vertices but of course not necessary for a regular vertex)
*/
struct HybridVertex
{
    HybridVertex(){} 
    HybridVertex(Intersect& i);
    HybridVertex(Polygon::Vertex v);
    glm::vec2 point; // only a sort of cache - can be calculated any time
    union {
        Polygon* owner;
        Polygon* edge1_owner;
    };
    union {
        int vertex;
        int edge1_index;
    };
    union {
        float alpha;
        float alpha1;
    };

    // For intersects only:
    bool intersect;
    Polygon* edge2_owner;
    int edge2_index;
    float alpha2;
};


struct IntersectionContact {
    float depth;
    glm::vec2 normal;
    EdgePoint subj_point;
    EdgePoint ref_point;
};

struct Intersection
{
    std::vector<HybridVertex> vertices;

    Intersection() {};
    Intersection(Polygon& p); // 'copy' p

    float signed_area();
    glm::vec2 centroid();
    glm::vec2 get_point(int vertex_number, float alpha);

    /* Direction to move not_reference back from the intersection */
    glm::vec2 find_normal_wrt(Polygon* not_reference, int start_vertex, int end_vertex);
    glm::vec2 find_normal_wrt(Polygon* not_reference);

    IntersectionContact get_contact(Polygon* reference, bool ref_outside, Polygon* subject, bool subj_outside, Renderer& renderer);

    Polygon* edge_owner(int edge_start_index);

    LineStrip cast_internal_shadow(glm::vec2 direction, Polygon* subject, Renderer &renderer);

    

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
    Polygon * find_parent_of_intersection_edge(Vertex v1, Vertex v2);

};




/*** Templates ***/

template <Intersect::Which which> // TODO separate functions for FIRST and SECOND
bool Intersect::lt(Intersect& i, Intersect& j)
{
    if (which == FIRST) {
        if (i.edge1.get_index() < j.edge1.get_index()) {
            return true;
        } else if (i.edge1.get_index() == j.edge1.get_index()) {
            return (i.alpha1 < j.alpha1);
        } else {
            return false;
        }
    } else {
        if (i.edge2.get_index() < j.edge2.get_index()) {
            return true;
        } else if (i.edge2.get_index() == j.edge2.get_index()) {
            return (i.alpha2 < j.alpha2);
        } else {
            return false;
        }
    }
}
