#pragma once

#include "constants.h"
#include "glutils.h"
#include "EdgePoint.h"


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
    template <typename IteratorType>
    class RangeAccessor;
    template <bool transformed>
    class EdgeIterator;
    template <bool transformed>
    class EdgeRangeIterator;
    template <bool transformed>
    class VertexIterator;
    template <bool transformed>
    class VertexRangeIterator;
 /** METHODS **/
 public:
	Polygon();
    Polygon(std::vector<glm::vec2> vertices);

    /* Iterating: */
    // Bad code ahead, should probably be moved to c++ file (not in class)
    Accessor<EdgeIterator<true>>
        edges() { return Accessor<EdgeIterator<true>>(*this); }
    RangeAccessor<EdgeRangeIterator<true>>
        edges(EdgePoint& start, EdgePoint& end) { return edges(start.index, start.alpha, end.index, end.alpha); }
    RangeAccessor<EdgeRangeIterator<true>>
        edges(int start_index, float start_alpha, int end_index, float end_alpha) { return RangeAccessor<EdgeRangeIterator<true>>(*this, start_index, start_alpha, end_index, end_alpha); }

    Accessor<EdgeIterator<false>>
        model_edges() { return Accessor<EdgeIterator<false>>(*this); }
    RangeAccessor<EdgeRangeIterator<false>>
        model_edges(EdgePoint& start, EdgePoint& end) { return model_edges(start.index, start.alpha, end.index, end.alpha); }
    RangeAccessor<EdgeRangeIterator<false>>
        model_edges(int start_index, float start_alpha, int end_index, float end_alpha) { return RangeAccessor<EdgeRangeIterator<false>>(*this, start_index, start_alpha, end_index, end_alpha); }

    Accessor<VertexIterator<true>>
        vertices() { return Accessor<VertexIterator<true>>(*this); }
    RangeAccessor<VertexRangeIterator<true>>
        vertices(EdgePoint& start, EdgePoint& end) { return vertices(start.index, start.alpha, end.index, end.alpha); }
    RangeAccessor<VertexRangeIterator<true>>
        vertices(int start_index, float start_alpha, int end_index, float end_alpha) { return RangeAccessor<VertexRangeIterator<true>>(*this, start_index, start_alpha, end_index, end_alpha); }

    Accessor<VertexIterator<false>>
        model_vertices() { return Accessor<VertexIterator<false>>(*this); }
    RangeAccessor<VertexRangeIterator<false>>
        model_vertices(EdgePoint& start, EdgePoint& end) { return model_vertices(start.index, start.alpha, end.index, end.alpha); }
    RangeAccessor<VertexRangeIterator<false>>
        model_vertices(int start_index, float start_alpha, int end_index, float end_alpha) { return RangeAccessor<VertexRangeIterator<false>>(*this, start_index, start_alpha, end_index, end_alpha); }

    /* Old iterating */
    Edge first_edge();
    Edge last_edge();
	int num_edges();
	LineSegment get_edge(int index);


    /* Access & Transformations */
    glm::vec2 transform(glm::vec2 point); // Transform from model to world coordinates
    glm::vec2 detransform(glm::vec2 point); // Transform from world to model coordinates
    glm::vec2 transformed(int vertex_index);
    glm::vec2 transform_center(glm::vec2 point, glm::vec2 center);

    glm::vec2 model_vertex(int index) { return _vertices[index]; }
    glm::vec2 edge_point(int index, float alpha);
    glm::vec2 model_edge_point(int index, float alpha);

    /* Info */
    int num_vertices()              { return _vertices.size(); }
    float get_moment_of_inertia()   { return moment_of_inertia; }
    float get_mass()                { return mass; }
    glm::vec2 get_center_of_mass()  { return center_of_mass; }
    bool is_CCW()                   { return CCW; }
    float get_radius()              { return radius; }
    AABB<2> calc_bounding_box();
    /* Query */
    bool vertex_is_concave(int index);



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













    /////////////////////////
    //
    // Iteration
    //
    /////////////////////////
    

    template <bool transformed>
    class EdgeIterator
    {
     public:
        ::Edge operator* () { return ::Edge(edge, get_index()); }
        EdgeIterator<transformed> & operator++ ();
        // TODO don't know if the following should be added //
        /* EdgeIterator<transformed> & operator-- (); */
        /* void set_index(int index)               {  this->index = index; } */
        bool operator== (EdgeIterator& other)   { return index == other.index; }
        bool operator!= (EdgeIterator& other)   { return index != other.index; }
        int get_index() {
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
    class EdgeRangeIterator
    {
     public:
        ::Edge operator* () { return ::Edge(edge, get_index()); }
        EdgeRangeIterator<transformed> & operator++ ();
        // TODO don't know if the following should be added //
        /* EdgeIterator<transformed> & operator-- (); */
        /* void set_index(int index)               {  this->index = index; } */
        bool operator== (EdgeRangeIterator& other)   { return index == other.index; }
        bool operator!= (EdgeRangeIterator& other)   { return index != other.index; }
        int get_index() {
            return (index - 1 + polygon._vertices.size()) % polygon._vertices.size();
        }
     private:
        EdgeRangeIterator(int index, float alpha, int end_index, float end_alpha, Polygon& polygon);
        EdgeRangeIterator(Polygon& polygon) : polygon(polygon), index(polygon.num_vertices()) {}; // constructs the end
        Polygon& polygon;
        int   index; // index refers to the _end_ of the edge
        float alpha;
        int   end_index;
        float end_alpha;
        LineSegment edge;

     template <typename IteratorType>
     friend class RangeAccessor;
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


    template <bool transformed>
    class VertexRangeIterator
    {
     public:
        ::Vertex operator* () { return ::Vertex(point, index); }
        VertexRangeIterator<transformed> & operator++ ();
        bool operator== (VertexRangeIterator& other) { return index == other.index; }
        bool operator!= (VertexRangeIterator& other) { return index != other.index; }
     private:
        VertexRangeIterator(int index, float alpha, int end_index, float end_alpha, Polygon& polygon);
        VertexRangeIterator(Polygon& polygon) : polygon(polygon), index(polygon.num_vertices()) {}; // constructs the end
        Polygon& polygon;
        int   index;
        float alpha;
        int   end_index;
        float end_alpha;
        glm::vec2 point;

     template <typename IteratorType>
     friend class RangeAccessor;
    };


    template <typename IteratorType>
    class Accessor
    {
     public:
        //// Note: end index doesn't matter for non-range-based iterators ////
        Accessor(Polygon& polygon)
            : polygon(polygon) {};
        /* Init iteration */
        IteratorType begin() {
            return IteratorType(false, polygon);
        }
        IteratorType end() {
            return IteratorType(true, polygon);
        }
     private:
        Polygon& polygon;
    };


    template <typename IteratorType>
    class RangeAccessor
    {
     public:
        //// Note: end index doesn't matter for non-range-based iterators ////
        RangeAccessor(Polygon& polygon, int start_index, float start_alpha, int end_index, float end_alpha)
            : polygon(polygon), start_index(start_index), start_alpha(start_alpha), end_index(end_index), end_alpha(end_alpha) {};
        RangeAccessor(Polygon& polygon, EdgePoint& start_point, EdgePoint& end_point)
            : polygon(*start_point.parent), start_index(start_point.index), start_alpha(start_point.alpha), end_index(end_point.index), end_alpha(end_point.alpha) {};
        /* Init iteration */
        IteratorType begin()
        {
            return IteratorType(start_index, start_alpha, end_index, end_alpha, polygon);
        }
        IteratorType end()
        {
            return IteratorType(polygon);
        }
     private:
        Polygon& polygon;
        int     start_index = 0;
        float   start_alpha = 0;
        int     end_index = 0;
        float   end_alpha = 0;
    };
};




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

    template <bool _transformed>
inline Polygon::EdgeIterator<_transformed>::EdgeIterator(bool is_end, Polygon& polygon)
    : polygon(polygon)
{
    if (is_end) { // end iterator
        this->index = polygon.num_vertices();
    } else {
        this->index = 0;
        if ( !_transformed ) {
            edge.start = polygon._vertices[polygon.num_vertices() - 1];
            edge.end = polygon._vertices[0];
        } else {
            edge.start = polygon.transformed(polygon.num_vertices() - 1);
            edge.end = polygon.transformed(0);
        }
    }
}

    template <bool _transformed>
inline Polygon::EdgeIterator<_transformed>& Polygon::EdgeIterator<_transformed>::operator++ ()
{
    assert(index >= 0);
    assert(index < polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size()) {
        edge.start = edge.end;
        if ( !_transformed )
            edge.end = polygon._vertices[index];
        else
            edge.end = polygon.transformed(index);
    }
    return *this;
}

////////////////////////////////////////////////////
//
// Edge Range Iterator
//
/////////////////////////////////////////////////////

    template <bool _transformed>
inline Polygon::EdgeRangeIterator<_transformed>::EdgeRangeIterator(int index, float alpha, int end_index, float end_alpha, Polygon& polygon)
    :polygon(polygon), index(index), alpha(alpha), end_index(end_index), end_alpha(end_alpha)
{
    if ( index != polygon.num_vertices() ) { // end iterator
        if ( !_transformed )
            edge.end = polygon.model_edge_point(index, alpha);
        else
            edge.end = polygon.edge_point(index, alpha);
        operator++ ();
    }
}

    template <bool _transformed>
inline Polygon::EdgeRangeIterator<_transformed>& Polygon::EdgeRangeIterator<_transformed>::operator++ ()
{
    assert(index != polygon.num_vertices());
    if (index == end_index) {
        if (alpha == 0) {
            alpha = end_alpha;
        } else {
            index = polygon.num_vertices(); // set to end if it went from the end vertex
            return *this;
        }
    } else {
        index = (index + 1) % polygon.num_vertices();
    }

    edge.start = edge.end;

    if ( !_transformed )
        edge.end = polygon.model_edge_point(index, alpha);
    else
        edge.end = polygon.edge_point(index, alpha);

    if (index == end_index)
        alpha = 0;
    return *this;
}

////////////////////////////////////////////////////
//
// Vertex iterator
//
/////////////////////////////////////////////////////

    template<bool _transformed>
inline Polygon::VertexIterator<_transformed>::VertexIterator(bool is_end, Polygon& polygon)
    :polygon(polygon)
{
    if (is_end)
        this->index = polygon.num_vertices();
    else {
        this->index = 0;
        if ( !_transformed )
            point = polygon._vertices[0];
        else
            point = polygon.transformed(0);
    }
}
    template<bool _transformed>
inline Polygon::VertexIterator<_transformed>& Polygon::VertexIterator<_transformed>::operator++ ()
{
    assert(index != polygon._vertices.size());
    ++ index;
    if (index < polygon._vertices.size()) {
        if ( !_transformed )
            point = polygon._vertices[index];
        else
            point = polygon.transformed(index);
    }
    return *this;
}

////////////////////////////////////////////////////
//
// Vertex Range Iterator
//
/////////////////////////////////////////////////////

template<bool _transformed>
inline Polygon::VertexRangeIterator<_transformed>::VertexRangeIterator(int index, float alpha, int end_index, float end_alpha, Polygon& polygon)
    :polygon(polygon), index(index), alpha(alpha), end_index(end_index), end_alpha(end_alpha)
{
    if ( index != polygon.num_vertices() ) {
        if ( !_transformed )
            point = polygon.model_edge_point(index, alpha);
        else
            point = polygon.edge_point(index, alpha);
    }
}
template<bool _transformed>
inline Polygon::VertexRangeIterator<_transformed>& Polygon::VertexRangeIterator<_transformed>::operator++ ()
{
    assert(index != polygon.num_vertices());
    std::cout << "Alpha: " << alpha << std::endl;
    if (index == end_index) {
        if (alpha == 0) {
            alpha = end_alpha;
        } else {
            index = polygon.num_vertices(); // set to end if it went from the end vertex
            return *this;
        }
    } else {
        index = (index + 1) % polygon.num_vertices();
    }
    if ( !_transformed )
        point = polygon.model_edge_point(index, alpha);
    else
        point = polygon.edge_point(index, alpha);

    if (index == end_index)
        alpha = 0;

    return *this;
}
