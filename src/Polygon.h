#pragma once

// class POLYGON

#include "glutils.h"
#include <glm/glm.hpp>
#include <vector>
#include <iterator>
#include <utility>
#include <iostream>
class SubPolygon;

struct Intersect;

struct NewVertex;
struct FullIntersect;

struct HybridVertex;
struct Intersection;

enum Axis { X=0, Y=1 };

typedef bool Direction;
const bool BACK = false;
const bool FORTH = true;

typedef std::pair<glm::vec2, glm::vec2> LineSegment;
struct Triangle
{
    Triangle(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec3 color):a(a), b(b), c(c), color(color) {};
    glm::vec2 a,b,c;
    glm::vec3 color;
};

class Polygon
{
public:
    static std::vector<Intersect> overlaps(Polygon& a, Polygon& b);
    static std::vector<Intersection> ExtractIntersections(Polygon& p, Polygon& q, bool flip_logic);
public:
	Polygon();

	std::vector<glm::vec2> vertices;
    glm::mat3 matrix;
    glm::vec2 transform(glm::vec2 point); // Transform from model to world coordinates
    glm::vec2 transformed(int vertex_index);
    glm::vec2 transform_center(glm::vec2 point, glm::vec2 center);


    void appendStencilTriangles(BufferWriter<float> &buffer);
    void appendLinesToVector(std::vector<float> &list);
	// SHAPE
	float signedArea();
	glm::vec2 centroid();
	float radius();
    glm::vec2 getPoint(int vertex_number, float alpha);

	// Monotonize, Triangulate, decompose into convex pieces --- returns #diagonals that are from step 1
	int decompose(std::vector<Triangle> &triangles, std::vector<LineSegment> &addedLines);




	// Iterating
	int numEdges();
	LineSegment getEdge(int index);

    ///////////////////////
    // Vertex smart-pointer
    ///////////////////////
    class Vertex
    {
    public:
        Vertex(): index(0), parent(0) {};
        Vertex(int index, Polygon *parent);
        int getIndex() { return index; }
        void setIndex(int val);
        Polygon* getParent() { return parent; }
        
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
        // bool operator< (Edge other);
        bool operator== (Edge other);
        int operator() (int x) const; // Gives y value at the given x value
        glm::vec2 & start() const;
        glm::vec2 & end() const;
        glm::vec2 start_tr() const; // World coordinates
        glm::vec2 end_tr() const;
        int getIndex() const { return index; }
        Polygon* getParent() const { return parent; }

    private:
        int index;
        Polygon *parent;
    };
    //////////////////////
    // Point on polygon //
    //////////////////////
    class Point
    {
    public:
        Point():start_index(0), alpha(0), parent(0) {};
        Point(int start_index, float alpha, Polygon *parent);
        bool operator== (Point other);
        glm::vec2 & start() const;
        glm::vec2 & end() const;
        glm::vec2 start_tr() const; // World coordinates
        glm::vec2 end_tr() const;
        int getIndex() const { return start_index; }
        Polygon* getParent() const { return parent; }
        float getAlpha() const { return alpha; }
    private:
        int start_index;
        float alpha; // how far from start to end the point is along the edge
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
        glm::vec2& start() const { return parent->vertices[start_index];}
        glm::vec2& end() const   { return parent->vertices[end_index];}
        int startIndex()        { return start_index;}
        int endIndex()          { return end_index;}

        Diagonal& operator= (Diagonal rhs) { parent=rhs.parent; start_index=rhs.start_index; end_index=rhs.end_index; return *this;}
    private:
        Polygon *parent;
        int start_index;
        int end_index;
    };

private:
	// Reverse: Go from right to left.
    // output: diagonals
	void monotonize(std::vector<Diagonal> &diagonals, bool reverse);
    // input: parts
    void triangulate(std::vector<SubPolygon> &parts, std::vector<Diagonal> &diagonals, std::vector<Triangle> &triangles);
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

	float signedArea();
    void fill(); // Let subpolygon be the full polygon
    bool containsDiagonal(int a, int b); // Returns true if vertices a and b are contained
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
        int getIndex() { return parent->indices[index]; }
        int getIndexIndex() { return index; }
        void setIndex(int val) { index = val; }
        
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
        int getIndex() const;
        int getIndexIndex() const { return index; }
        SubPolygon& getParent() const { return *parent; }
    private:
        int index;
        SubPolygon *parent;
    };
};
std::ostream& operator<< (std::ostream& lhs, SubPolygon& p);

