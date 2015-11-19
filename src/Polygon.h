#pragma once

// class POLYGON

#include "Geometry.h"

#include <glm/glm.hpp>
#include <vector>
#include <iterator>
#include <utility>

class SubPolygon;

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
	Polygon();

	glm::vec3 color;
	std::vector<glm::vec2> vertices;

	// SHAPE
	float signedArea();
	glm::vec2 centroid();
	float radius();

	// Monotonize, Triangulate, decompose into convex pieces
	std::vector<LineSegment> decompose(std::vector<Triangle> &triangles);


	// Iterating
	int numEdges();
	LineSegment getEdge(int index);

    ///////////////////////
    // Vertex smart-pointer
    ///////////////////////
    class Vertex
    {
    public:
        Vertex(int index, Polygon &parent);
        int getIndex() { return index; }
        void setIndex(int val) { index = val; }
        
        glm::vec2 & operator* ();
        glm::vec2 * operator-> ();
        glm::vec2 & preceding();
        glm::vec2 & successive();

        Vertex& operator= (Vertex rhs) { parent=rhs.parent; index=rhs.index; return *this;}


    private:
        int index;
        Polygon &parent;
    };

    ///////////////////////
    // Edge smart-pointer
    ///////////////////////
    class Edge
    {
    public:
        Edge(int index, Polygon &parent);
        // bool operator< (Edge other);
        bool operator== (Edge other);
        int operator() (int x) const; // Gives y value at the given x value
        glm::vec2 & start() const;
        glm::vec2 & end() const;
        int getIndex() const { return index; }

        Edge& operator= (Edge rhs) { parent=rhs.parent; index=rhs.index; return *this;}
    private:
        int index;
        Polygon &parent;
    };

    /////////////////////////
    // Diagonal smart-pointer
    /////////////////////////
    class Diagonal
    {
    public:
        Diagonal(int vertex1, int vertex2, Polygon &parent)
            :start_index(vertex1), end_index(vertex2), parent(parent) {};
        glm::vec2& start() const { return parent.vertices[start_index];}
        glm::vec2& end() const   { return parent.vertices[end_index];}
        int startIndex()        { return start_index;}
        int endIndex()          { return end_index;}

        Diagonal& operator= (Diagonal rhs) { parent=rhs.parent; start_index=rhs.start_index; end_index=rhs.end_index; return *this;}
    private:
        Polygon &parent;
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

bool overlaps(Polygon a, Polygon b);


class SubPolygon
{
public:

    SubPolygon(Polygon &mother):mother(mother) {};
    std::vector<int> indices; // Indices in polygon vertex buffer
    Polygon& mother;

    void fill(); // Let subpolygon be the full polygon
    bool containsDiagonal(int a, int b); // Returns true if vertices a and b are contained
    void split(int a, int b, SubPolygon& out1, SubPolygon& out2); // Split into two polygons at edge (a, b)
    // Rule of three
    SubPolygon& operator= (SubPolygon p);

    ///////////////////////
    // Vertex smart-pointer
    ///////////////////////
    class Vertex
    {
    public:
        // index: This is the index into the indices array of the subpolygon
        Vertex(int index, SubPolygon &parent);
        int getIndex() { return parent.indices[index]; }
        int getIndexIndex() { return index; }
        void setIndex(int val) { index = val; }
        
        glm::vec2 & operator* ();
        glm::vec2 * operator-> ();
        glm::vec2 & preceding();
        glm::vec2 & successive();

        Vertex& operator= (Vertex rhs) { parent=rhs.parent; index=rhs.index; return *this;}
    private:
        int index;
        SubPolygon &parent;
    };
};
std::ostream& operator<< (std::ostream& lhs, SubPolygon& p);

