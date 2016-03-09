#pragma once
#include "Polygon.h"

struct Intersect;

struct NewVertex;
struct FullIntersect;

struct HybridVertex;
struct Intersection;

// External:
class Polygon;

typedef bool Side;
const bool IN = true; // = FORTH
const bool OUT = false; // = BACK


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


struct Manifold {
    glm::vec2 normal; // includes depth

    Polygon* reference;
    glm::vec2 ref_point;
    Polygon* subject;
    glm::vec2 subj_point;

};

struct Intersection // Or HybridPolygon..
{
    std::vector<HybridVertex> vertices;

    Intersection() {};
    Intersection(Polygon& p); // 'copy' p
    Intersection ExtendInDirection(glm::vec2 v, Polygon* subject);
    Manifold manifold(glm::vec2 relative_velocity, Polygon* subject, Polygon* reference);
    float signedArea();
    glm::vec2 centroid();
    void appendLinesToVector(std::vector<float> &list);
    glm::vec2 getPoint(int vertex_number, float alpha);

    class Vertex {
    public:
        Vertex(): index(0), parent(0) {};
        Vertex(int index, Intersection *parent);
        int getIndex() { return index; }
        void setIndex(int val);
        Intersection* getParent() { return parent; }
        
        HybridVertex & operator* ();
        HybridVertex * operator-> ();

        // Changes what vertex is pointed to
        Vertex& operator++ ();
        Vertex& operator-- ();
        //
        bool operator== (Vertex& v);

        HybridVertex & preceding();
        HybridVertex & successive();

    private:
        int index;
        Intersection *parent;
    };
};


struct NewVertex {
    NewVertex(Intersect& i, Side in_out) {
        point = i.point;
        intersect = true;
        processed = false;
        this->in_out = in_out;

        vertex = HybridVertex(i);
    }
    NewVertex(Polygon::Vertex vert, Side in_out) {
        point = vert.transformed();
        intersect = false;
        this->in_out = in_out;

        vertex = HybridVertex(vert);
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

struct FullIntersect { // An intersection is needed to sort the NewVertices
    FullIntersect(Intersect *i, NewVertex *p, NewVertex *q): i(i), vert_p(p), vert_q(q){}
    Intersect* i;
    NewVertex* vert_p;
    NewVertex* vert_q;
    template <Intersect::Which which> // sort wrt p or q?
    static bool lt(FullIntersect& a, FullIntersect& b);
};


