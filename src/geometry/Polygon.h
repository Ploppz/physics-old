#pragma once

#include "constants.h"
#include "glutils.h"


#include <glm/glm.hpp>
#include <vector>
#include <iterator>
#include <utility>
#include <iostream>
#include <tuple>

class SubPolygon;

struct Intersect;

struct HybridVertex;
struct Intersection;

template <int d>
struct AABB;


struct LineSegment {
    LineSegment() : start{}, end{} {}
    LineSegment(glm::vec2 start, glm::vec2 end) : start(start), end(end) {}
    glm::vec2 start;
    glm::vec2 end;
};

/* Used for iterating:
 * I don't like it so much, waiting for C++17 structured bindings
 * in order to return std::tuple<LineSegment&, int>. */
struct Edge {
    Edge(LineSegment& edge,int index) : start(edge.start), end(edge.end), index(index) {}
    glm::vec2 start;
    glm::vec2 end;
    int index;
};
struct Vertex {
    Vertex(glm::vec2 point, int index) : point(point), index(index) {}
    glm::vec2 point;
    int index;
};

struct Triangle
{
    Triangle(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec3 color):a(a), b(b), c(c), color(color) {};
    glm::vec2 a,b,c;
    glm::vec3 color;
};

/* TODO not used yet, but could be useful? */
struct Transform
{
    Transform(glm::vec2 translation, float orientation) : translation(translation), orientation(orientation) {}
    glm::vec2 translation;
    float orientation;
};

class Polygon
{
 /** STATIC **/
 public:
    static std::vector<Intersect> find_intersects(Polygon& a, Polygon& b);
    static std::vector<Intersection> extract_intersections(Polygon& p, Polygon& q, bool flip_p_logic, bool flip_q_logic);

    friend SubPolygon;
    class Edge;
    class Diagonal;
    template <typename IteratorType>
    class Accessor;
    template <bool transformed>
    class EdgeIterator;
    template <bool transformed>
    class VertexIterator;
 /** METHODS **/
 public:
	Polygon();
    Polygon(std::vector<glm::vec2> vertices);

    /* Iterating: */
    Accessor<EdgeIterator<true>> edges() { return Accessor<EdgeIterator<true>>(*this); }
    Accessor<EdgeIterator<false>> model_edges() { return Accessor<EdgeIterator<false>>(*this); }
    Accessor<VertexIterator<true>> vertices() { return Accessor<VertexIterator<true>>(*this); }
    Accessor<VertexIterator<false>> model_vertices() { return Accessor<VertexIterator<false>>(*this); }
    /* Old iterating */
    Edge first_edge();
    Edge last_edge();
	int num_edges();
	LineSegment get_edge(int index);


    /* Access & Transformations */
    glm::vec2 transform(glm::vec2 point); // Transform from model to world coordinates
    glm::vec2 detransform(glm::vec2 point); // Transform from world to model coordinates
    glm::vec2 transformed(int vertex_index);
    glm::vec2 model_vertex(int index) { return _vertices[index]; }
    glm::vec2 transform_center(glm::vec2 point, glm::vec2 center);

    /* Info */
    int num_vertices()              { return _vertices.size(); }
    float get_moment_of_inertia()   { return moment_of_inertia; }
    float get_mass()                { return mass; }
    glm::vec2 get_center_of_mass()  { return center_of_mass; }
    bool is_CCW()                   { return CCW; }
    float get_radius()              { return radius; }
    AABB<2> calc_bounding_box();



	// Monotonize, Triangulate, decompose into convex pieces --- returns #diagonals that are from step 1
	int decompose(std::vector<Triangle> &triangles, std::vector<LineSegment> &added_lines);
 private:
	void monotonize(std::vector<Diagonal> &out_diagonals, bool reverse);
    void triangulate(std::vector<SubPolygon> &parts, std::vector<Diagonal> &diagonals, std::vector<Triangle> &triangles);

    float calculate_moment_of_inertia();
	float signed_area();
    void apply_center_of_mass(glm::vec2 center);
	float calculate_radius();

    void calculate_shape_dependent_variables();

 /** MEMBERS **/
 public:
    glm::vec2 position;
    float orientation;
    
 private:
	std::vector<glm::vec2> _vertices;
    /* Shape-dependent calculations: */
    float moment_of_inertia;
    float mass;
    float radius;
    glm::vec2 center_of_mass;
    bool CCW;

public: /** Helper classes **/

    ///////////////////////
    // Vertex smart-pointer
    ///////////////////////
    class Vertex
    {
    public:
        Vertex(): index(0), parent(0) {};
        Vertex(int index, Polygon *parent);
        int get_index() { return index; }
        void set_index(int val);
        Polygon* get_parent() { return parent; }
        
        glm::vec2 & operator* ();
        glm::vec2 * operator-> ();

        // Changes what vertex is pointed to
        Vertex& operator++ ();
        Vertex& operator-- ();
        bool operator== (Vertex& v);

        glm::vec2 & preceding();
        glm::vec2 & successive();

        // Vertex transformed WRT to parent polygon
        glm::vec2 transformed();

    private:
        int index;
        Polygon *parent;
    };

    ///////////////////////
    // Edge smart-pointer
    ///////////////////////
    class Edge
    {
    public:
        Edge():index(0), parent(0) {};
        Edge(int index, Polygon *parent);

        glm::vec2 normal_tr();

        int operator() (int x) const; // Gives y value at the given x value
        glm::vec2 & start() const;
        glm::vec2 & end() const;
        glm::vec2 start_tr() const; // World coordinates
        glm::vec2 end_tr() const;
        int get_index() const { return index; }
        Polygon* get_parent() const { return parent; }

        bool operator== (Edge other); // Asumes that they have the same polygon!
        bool operator!= (Edge other);
        // Changes what edge is pointed to
        Edge& operator++ ();
        Edge& operator-- ();

    private:
        int index;
        Polygon *parent;
    };

    template <bool transformed>
    class EdgeIterator
    {
     public:
        ::Edge operator* () { return ::Edge(edge, get_index()); }
        EdgeIterator<transformed> & operator++ ();
        bool operator== (EdgeIterator& other) { return index == other.index; }
        bool operator!= (EdgeIterator& other) { return index != other.index; }
        int get_index()
        {
            return (index - 1 + polygon._vertices.size()) % polygon._vertices.size();
        }
     private:
        EdgeIterator(bool is_end, Polygon& polygon);
        Polygon& polygon;
        int index; // index refers to the _end_ of the edge
        LineSegment edge;

    template <typename IteratorType>
     friend class Accessor;

    };

    template <bool transformed>
    class VertexIterator
    {
     public:
        ::Vertex operator* () { return ::Vertex(point, index); }
        VertexIterator<transformed> & operator++ ();
        bool operator== (VertexIterator& other) { return index == other.index; }
        bool operator!= (VertexIterator& other) { return index != other.index; }
     private:
        VertexIterator(bool is_end, Polygon& polygon);
        Polygon& polygon;
        int index;
        glm::vec2 point;

    template <typename IteratorType>
     friend class Accessor;

    };

    template <typename IteratorType>
    class Accessor
    {
     public:
        Accessor(Polygon& polygon) : polygon(polygon) {};
        /* Init iteration */
        IteratorType begin()
        {
            return IteratorType(false, polygon);
        }
        IteratorType end()
        {
            return IteratorType(true, polygon);
        }
     private:
        Polygon& polygon;
    };

    /////////////////////////
    // Diagonal smart-pointer
    /////////////////////////
    class Diagonal
    {
    public:
        Diagonal(int vertex1, int vertex2, Polygon *parent)
            :parent(parent), start_index(vertex1), end_index(vertex2) {};
        glm::vec2& start() const { return parent->_vertices[start_index];}
        glm::vec2& end() const   { return parent->_vertices[end_index];}

        Diagonal& operator= (Diagonal rhs) { parent=rhs.parent; start_index=rhs.start_index; end_index=rhs.end_index; return *this;}
// was private ..
        Polygon *parent;
        int start_index;
        int end_index;
    };
};


/*********************/
/*    SubPolygon     */
/*********************/


class SubPolygon
{
public:
    SubPolygon():mother(0) {};
    SubPolygon(Polygon *mother):mother(mother) {};
    std::vector<int> indices; // Indices in polygon vertex buffer
    Polygon* mother;

	float signed_area();
    void fill(); // Let subpolygon be the full polygon
    bool contains_diagonal(int a, int b); // Returns true if vertices a and b are contained
    void split(int a, int b, SubPolygon& out1, SubPolygon& out2); // Split into two polygons at edge (a, b)
    template <Axis axis>
    void monotonize(std::vector<Polygon::Diagonal> &diagonals, Direction dir);

    ///////////////////////
    // Vertex smart-pointer
    ///////////////////////
    class Vertex
    {
    public:
        // index: This is the index into the indices array of the subpolygon
        Vertex() :index(0), parent(0) {}
        Vertex(int index, SubPolygon *parent);
        int get_index() { return parent->indices[index]; }
        int get_index_index() { return index; }
        void set_index(int val) { index = val; }
        
        glm::vec2 & operator* ();
        glm::vec2 * operator-> ();
        glm::vec2 & preceding();
        glm::vec2 & successive();
    private:
        int index;
        SubPolygon *parent;
    };

    ///////////////////////
    // Edge smart-pointer
    ///////////////////////
    class Edge
    {
    public:
        Edge() :index(0), parent(0) {}
        Edge(int index, SubPolygon *parent);
        // bool operator< (Edge other);
        bool operator== (Edge other) const;
        int y(int x) const;
        int x(int y) const;
        glm::vec2 & start() const;
        glm::vec2 & end() const;
        int get_index() const;
        int get_index_index() const { return index; }
        SubPolygon& get_parent() const { return *parent; }
    private:
        int index;
        SubPolygon *parent;
    };
};
std::ostream& operator<< (std::ostream& lhs, SubPolygon& p);




////////////////////////////////////////////////////
// Thoughts on iterators
//  - Accessor(start_index, end_index). The iterator makes sure that after the end_index, index = some sentinel
//  - With the current system, index always points to the actual index, but we have to do a test after incrementing
//      - Would it be better to just start at index=-1, and get_index() = index + 1?
/////////////////////////////////////////////////////
//
// Edge iterator
//
/////////////////////////////////////////////////////

template <>
inline Polygon::EdgeIterator<false>::EdgeIterator(bool is_end, Polygon& polygon)
    : polygon(polygon)
{
    if (is_end) { // end iterator
        this->index = polygon._vertices.size();
    } else {
        this->index = 0;
        edge.start = polygon._vertices[polygon._vertices.size() - 1];
        edge.end = polygon._vertices[0];
    }
}
template <>
inline Polygon::EdgeIterator<true>::EdgeIterator(bool is_end, Polygon& polygon)
    : polygon(polygon)
{
    if (is_end) { // end iterator
        this->index = polygon._vertices.size();
    } else {
        this->index = 0;
        edge.start = polygon.transformed(polygon._vertices.size() - 1);
        edge.end = polygon.transformed(0);
    }
}


template <>
inline Polygon::EdgeIterator<false>& Polygon::EdgeIterator<false>::operator++ ()
{
    assert(index >= 0);
    assert(index < polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size()) {
        edge.start = edge.end;
        edge.end = polygon._vertices[index];
    }
    return *this;
}
template <>
inline Polygon::EdgeIterator<true>& Polygon::EdgeIterator<true>::operator++ ()
{
    assert(index >= 0);
    assert(index < polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size()) {
        edge.start = edge.end;
        edge.end = polygon.transformed(index);
    }
    return *this;
}


////////////////////////////////////////////////////
//
// Vertex iterator
//
/////////////////////////////////////////////////////

template<>
inline Polygon::VertexIterator<false>::VertexIterator(bool is_end, Polygon& polygon)
    :polygon(polygon)
{
    if (is_end)
        this->index = polygon._vertices.size();
    else {
        this->index = 0;
        point = polygon._vertices[0];
    }
}
template<>
inline Polygon::VertexIterator<true>::VertexIterator(bool is_end, Polygon& polygon)
    :polygon(polygon)
{
    if (is_end)
        this->index = polygon._vertices.size();
    else {
        this->index = 0;
        point = polygon.transformed(0);
    }
}
template<>
inline Polygon::VertexIterator<false>& Polygon::VertexIterator<false>::operator++ ()
{
    assert(index != polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size())
        point = polygon._vertices[index];
    return *this;
}
template<>
inline Polygon::VertexIterator<true>& Polygon::VertexIterator<true>::operator++ ()
{
    assert(index != polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size())
        point = polygon.transformed(index);
    return *this;
}
