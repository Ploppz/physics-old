#include "Polygon.h"
#include "Geometry.h"
#include "tmp.h"

#include <list>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <set>

int abs_mod(int n, int range)
{
    return (n + range) % range;
}

Polygon::Polygon()
	: color(0.6f, 0.6f, 0.6f)
{
}


float Polygon::signedArea()
{
	float A = 0;
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end()) break;
		A += it->x * next->y - next->x * it->y;
	}
	A *= 0.5f;
	return A;
}
float SubPolygon::signedArea()
{
    float A = 0;
	std::vector<int>::iterator next;
	for (auto it = indices.begin(); it != indices.end(); it ++)
	{
		next = it; next ++;
		if (next == indices.end()) break;
		A += mother->vertices[*it].x * mother->vertices[*next].y - mother->vertices[*next].x * mother->vertices[*it].y;
	}
	A *= 0.5f;
	return A;
}

// LOCAL
glm::vec2 Polygon::centroid()
{
	float A = signedArea();
	glm::vec2 C {};
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end()) break;
		// Optimize this maybe..
		C.x += (it->x + next->x) * (it->x * next->y - next->x * it->y);
		C.y += (it->y + next->y) * (it->x * next->y - next->x * it->y);
	}
	C /= 6 * A;

	return C;
}
float Polygon::radius()
{
	glm::vec2 c = centroid();
	// Find point whose distance is MAXIMUM from c.
	float m = 0;
	float d;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		d = distance(c, *it); // TODO: made this relative to local coor system
		if (d > m) m = d;
	}
	return m;
}


// http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
// (3rd reply)


/* y-axis less than for Edges */
// 	We assume that the edges never cross, and that they do share some x-coordinate.
// 		Thus, we only need to find each of their left-most vertices, and see which one is the most right-most,
// 		and compare at this point.
// 	TODO: WARNING: does not tolerate vertical lines YET
template <Axis axis>
struct EdgeComparator
{
    bool operator() (const SubPolygon::Edge &lhs, const SubPolygon::Edge &rhs);
};
template <>
bool EdgeComparator<X>::operator() (const SubPolygon::Edge &lhs, const SubPolygon::Edge &rhs)
{
    if (lhs == rhs) return false;
	glm::vec2 subject, a_leftmost, b_leftmost;
	if (lhs.start().x < lhs.end().x) 	a_leftmost = lhs.start();
	else 						    	a_leftmost = lhs.end();
	if (rhs.start().x < rhs.end().x)    	b_leftmost = rhs.start();
	else 						    	b_leftmost = rhs.end();
	if (a_leftmost.x < b_leftmost.x) {
        // Do y-comparision at b's left-most vertex
		return lhs.y(b_leftmost.x) < b_leftmost.y;
	} else {
        // Do y-comparision at a's left-most vertex
		return a_leftmost.y < rhs.y(a_leftmost.x);
	}
}
template <>
bool EdgeComparator<Y>::operator() (const SubPolygon::Edge &lhs, const SubPolygon::Edge &rhs)
{
    if (lhs == rhs) return false;
	glm::vec2 subject, a_uppermost, b_uppermost;
	if (lhs.start().y < lhs.end().y) 	a_uppermost = lhs.start();
	else 						    	a_uppermost = lhs.end();
	if (rhs.start().y < rhs.end().y)    b_uppermost = rhs.start();
	else 						    	b_uppermost = rhs.end();
	if (a_uppermost.y < b_uppermost.y) {
        // Do y-comparision at b's left-most vertex
		return lhs.x(b_uppermost.y) < b_uppermost.x;
	} else {
        // Do y-comparision at a's left-most vertex
		return a_uppermost.x < rhs.x(a_uppermost.y);
	}
}

template <typename T>
void erase(std::set<SubPolygon::Edge, T>& set, int edge_index, SubPolygon* sub_polygon)
{
    auto it = set.find( SubPolygon::Edge(edge_index, sub_polygon));
    if (it != set.end()) set.erase(it);

}


template<Axis> struct axis_tag { };

//////////////////////////////
// Get Edge Just Below
/////////////////////////////
// With helpers to help partially instantiate the function templates

template <typename T_set>
SubPolygon::Edge getEdgeJustBelow(T_set s, SubPolygon::Vertex subject, axis_tag<X>)
{
    // X-monotonity: below on y-axis
	// Until later... linear search.
    auto next = s.begin();
	for (auto it = s.begin(); it != s.end(); it ++)
	{
		next = it; next ++;
		if (next == s.end()) next = s.begin();
		if ((*next).y(subject->x) > subject->y) {
			if ((*it).y(subject->x) > subject->y)
				std::cout << "getEdgeJustBelow(): Error: No available edges." << std::endl;
			else
				return *it;
		}
	}
	// Return last element
	assert(s.size() > 0);
    SubPolygon::Edge result = *s.rbegin();
	return result;
}
template <typename T_set>
SubPolygon::Edge getEdgeJustBelow(T_set s, SubPolygon::Vertex subject, axis_tag<Y>)
{
    // Y-monotonity: below on x-axis
	// Until later... linear search.
    auto next = s.begin();
	for (auto it = s.begin(); it != s.end(); it ++)
	{
		next = it; next ++;
		if (next == s.end()) next = s.begin();
		if ((*next).x(subject->y) > subject->x) {
			if ((*it).x(subject->y) > subject->x)
				std::cout << "getEdgeJustBelow(): Error: No available edges." << std::endl;
			else
				return *it;
		}
	}
	// Return last element
	assert(s.size() > 0);
    SubPolygon::Edge result = *s.rbegin();
	return result;
}

template <Axis axis, typename T_set>
SubPolygon::Edge getEdgeJustBelow(T_set s, SubPolygon::Vertex subject)
{
    return getEdgeJustBelow(s, subject, axis_tag<axis>());
}

////////////////////////
// Vertex comparision //
////////////////////////
//

template <typename T>
bool real_less_than(T a, T b, axis_tag<X>){return a->x < b->x; }
template <typename T>
bool real_less_than(T a, T b, axis_tag<Y>){ return a->y < b->y; }
template <Axis axis, typename T>
bool less_than(T a, T b){
    return real_less_than(a, b, axis_tag<axis>());
}


template <typename T>
bool real_greater_than(T a, T b, axis_tag<X>){ return a->x > b->x; }
template <typename T>
bool real_greater_than(T a, T b, axis_tag<Y>){ return a->y > b->y; }
template <Axis axis, typename T>
bool greater_than(T a, T b){
    return real_greater_than(a, b, axis_tag<axis>());
}
/////////////////////



//TODO WARNING: using a set for intersecting edges, edges with same y-value of start vertex is not allowed!

// Sorting: could sort a list of indices but complicates the comparision function

std::ostream &operator << (std::ostream &lhs, glm::vec2 &rhs);


int Polygon::decompose(std::vector<Triangle> &triangles, std::vector<LineSegment> &addedLines)
{
	std::vector<Diagonal> diagonals;

    // Make a full subpolygon which you make x-monotone, split into subpolygons which you make y-monotone
    SubPolygon full = SubPolygon(this);
    full.fill();
    full.monotonize<X>(diagonals, FORTH);
    full.monotonize<X>(diagonals, BACK);
    int split = diagonals.size();

    // SPLIT the polygon into subpolygons given the diagonals
        std::vector<SubPolygon> parts;
        SubPolygon entire(this); entire.fill();
        parts.push_back(entire);
        // Loop through diagonals and find out which sub polygons to further split
        SubPolygon p1(this), p2(this);
        for (auto d = diagonals.begin(); d != diagonals.end(); d ++)
        {
            for (int i = parts.size() - 1; i >= 0; i --)
            {
                if (parts[i].containsDiagonal(d->startIndex(), d->endIndex())) {
                    parts[i].split(d->startIndex(), d->endIndex(), p1, p2);
                    parts[i] = p1;
                    parts.push_back(p2);
                }
            }
        }
        auto part = parts.begin();
        {
            std::cout << *part << std::endl;
            triangulate(parts, diagonals, triangles);
        }

        std::cout << "\n\n\n" << std::endl;

    // Lines for drawing, pretty much
    for (auto it = diagonals.begin(); it != diagonals.end(); it ++)
    {
        // addedLines.push_back(LineSegment(it->start(), it->end()));
    }
    return split;
}

// TODO: Make static branches (template + structs holding different function?) instead of axis&reverse
// http://stackoverflow.com/questions/3615439/template-parameters-define-and-code-duplication
template <Axis axis>
void SubPolygon::monotonize(std::vector<Polygon::Diagonal> &diagonals, Direction dir)
{
    // std::cout << "Size: " << indices.size() << std::endl;
	bool CCW = (signedArea() > 0); // Counter clockwise
    bool logic = CCW ^ dir ^ axis; // CCW and dir flips the logic sometimes
	std::set<Edge, EdgeComparator<axis>> status;  // Updated sorted list of edges intersecting the scanline
        // (Only edges that have the polygon inside above)
	std::vector<Vertex> helpers(indices.size()); // A helper is the right-most vertex that an edge
		// can connect to with a vertical line
	std::vector<Vertex> events; // Sorted vertices by x-axis
	{// Sort vertices for scanline (keep in separate array)
		for (int i = 0; i < indices.size(); i ++)
		{
            // If the direction is backwards, vertices with higher index should preced vertices with lower index
            // If not flexible enough: sort by ascending or descending y-axis dependent on dir.
            if (dir == FORTH) {
                events.push_back(Vertex(i, this));
            } else {
                events.push_back(Vertex(indices.size() - 1 - i, this));
            }
		}
        auto fn = (dir == FORTH)? less_than<axis, Vertex> : greater_than<axis, Vertex>;
		std::sort(events.begin(), events.end(), fn);
	}
	
	/*
		A vertex can be a: start, end, split, merge or normal vertex
   */
	int index;
	std::set<Edge>::iterator it;
	float tmp;

	Direction precSide, succSide; // which side of the scanline the neighbouring points (preceding, succesive) are on
	bool concave; // is the current vertex concave or convex?
	glm::vec2 precRotated; // preceding vector, rotated 90 degrees.
	float dotProduct;
    std::cout << "EVENTS" << std::endl;
    for (auto it = events.begin(); it != events.end(); it ++) {
        std::cout << it->getIndex() << ", ";
    }
    std::cout << std::endl;

	for (int i = 0; i < events.size() - 1; i ++) // Each 'event[i]' is a vertex.
	{
        std::cout << "Vertex " << events[i].getIndex() << " .. " << events[i]->x << ", " << events[i]->y << std::endl;

		// True means that the preceding/successive vertex is BEHIND the scanline (Direction.BACK)
        if (axis == X) {
            precSide = Direction(events[i].preceding().x < events[i]->x);
            succSide = Direction(events[i].successive().x < events[i]->x);
            // Check for vertical edges
            if (events[i].preceding().x - events[i]->x == 0) precSide = succSide;
            if (events[i].successive().x - events[i]->x == 0) succSide = Direction(!bool(precSide));
        } else {
            precSide = Direction(events[i].preceding().y < events[i]->y);
            succSide = Direction(events[i].successive().y < events[i]->y);
            // Check for vertical edges
            if (events[i].preceding().y - events[i]->y == 0) precSide = succSide;
            if (events[i].successive().y - events[i]->y == 0) succSide = Direction(!bool(precSide));
        }

		// Used only if both vertices are on the same side of the scanline..:
		precRotated = events[i].preceding() - *events[i];
		tmp = precRotated.x;
		precRotated.x = - precRotated.y;
		precRotated.y = tmp;
		dotProduct = glm::dot(precRotated, events[i].successive() - *events[i]);
		concave  = (dotProduct < 0) ^ CCW;
		
        std::cout << std::boolalpha << "\t(" << bool(precSide) << ", " << bool(succSide) << std::endl;
        std::cout << "\t";
        if (precSide == succSide  &&  dir != precSide) // Both adjacent vertices are BEHIND the scanline
		{
			if (concave) { /* MERGE POINT */
                std::cout << "merge" << std::endl;
				// Remove the top-most edge from status
				index = events[i].getIndex() - (logic);
				if (index == -1) index = indices.size() - 1;
                erase(status, index, this);
				// Update helper of the edge that is just below this vertex
                Edge justBelow = getEdgeJustBelow<axis>(status, events[i]);
				helpers[justBelow.getIndex()] = events[i];
			} else { /* END POINT */
                std::cout << "end" << std::endl;
				// Remove the bottom-most edge from status
				index = events[i].getIndex() - (logic);
				if (index == -1) index = indices.size() - 1;
                erase(status, index, this);
			}
		} else if (precSide == succSide  &&  dir == precSide) // Both adjacent vertices are IN FRON OF the scanline
		{
			if (concave) { /* SPLIT POINT */
                std::cout << "split" << std::endl;
				// Find the first edge on the line projected down from this vertex
                Edge justBelow = getEdgeJustBelow<axis>(status, events[i]);
				// Make a diagonal to that edge's helper
                diagonals.push_back(
                        Polygon::Diagonal(events[i].getIndex(), helpers[justBelow.getIndex()].getIndex(), mother)  );
				// Update the helper of that edge
				helpers[justBelow.getIndex()] = events[i];
				// Insert the top-most adjacent edge to 'status', with this vertex as helper
				index = events[i].getIndex() - !(logic);
				if (index == -1) index = indices.size() - 1;
				status.insert( Edge(index, this) );
				helpers[index] = events[i];
			} else { /* START POINT */
                std::cout << "start" << std::endl;
				// Insert the bottom-most incident edge to 'status'
				index = events[i].getIndex() - !(logic);
				if (index == -1) index = indices.size() - 1;
                Edge e = Edge(index, this);
				status.insert( e);
				helpers[index] = events[i]; // with this vertex as helper
			}
		} else { /* NORMAL VERTEX */
            std::cout << "normal" << std::endl;
			// CW: If the incident edge (any of the two) goes to the right, then this vertex is in the _ceiling_
			bool ceiling;
            ceiling = (precSide == BACK && succSide == FORTH) ^ CCW;
            if (axis == X) {
                // ceiling = !(events[i]->x > events[i].successive().x) ^ CCW;
            } else {
                // ceiling = !(events[i]->y < events[i].successive().y) ^ CCW;
            }

			if (ceiling) {
                std::cout << "-> ceiling" << std::endl;
				// Replace helper of the edge just below
                Edge justBelow = getEdgeJustBelow<axis>(status, events[i]);
				helpers[justBelow.getIndex()] = events[i];
			} else {
                std::cout << "-> floor" << std::endl;
				// Remove edge frjm the status
				index = events[i].getIndex() - (logic);
				if (index == -1) index = indices.size() - 1;
                erase(status, index, this);
                // Add edge to status & update helper of the added edge
				index = events[i].getIndex() - !(logic);
				if (index == -1) index = indices.size() - 1;
				status.insert(Edge(index, this));
				helpers[index] = events[i];
			}
		}
        std::cout << "Status: ";
        for (auto it = status.begin(); it != status.end(); it ++)
        {
            std::cout << it->getIndex() << ", ";
        }
        std::cout << std::endl;
	}
    
}
void Polygon::triangulate(std::vector<SubPolygon> &parts, std::vector<Diagonal> &diagonals, std::vector<Triangle> &triangles)
{
	bool CCW = (signedArea() > 0); // Counter clockwise
    // Confusion: sub polygons have an indexed array of indices of the vertices
    // TRIANGULATE each subpolygon
    // We need to know which chain (upper or lower) each vertex is on. Indexed by index in original polygon.
    std::vector<bool> side;
    side.resize(vertices.size());
    //TODO Optimize - updating side info & sorting vertices could be done simultaneously
    for (auto part = parts.begin(); part != parts.end(); part ++) 
    {
    // Update chain (side) information
        std::vector<SubPolygon::Vertex> events;
        // Find min and max Vertices (x axis)
        Vertex max(part->indices[0], this), min(part->indices[0], this);
        int real_index, max_index=0, min_index=0; // Index in subpolygon
        for (unsigned int i = 1; i < part->indices.size(); i ++)
        {
            real_index = part->indices[i];
            if (vertices[real_index].x < min->x) {min.setIndex(real_index); min_index = i;}
            if (vertices[real_index].x > max->x) {max.setIndex(real_index); max_index = i;}
        }
        // Now that we know the subpolygon indices of the min and max vertices, we can split it up in two chains
        // Make sure that upper chain false and bottom chain true.
        bool current_side = CCW;
        for (int i = min_index; ; i ++) {
            if (i >= part->indices.size()) i = 0;
            if (i == max_index) break;
            side[part->indices[i]] = current_side;
            std::cout << "\tfalse " << part->indices[i] << std::endl;
        }
        current_side = !current_side;
        for (int i = max_index; ; i ++) {
            if (i >= part->indices.size()) i = 0;
            if (i == min_index) break;
            side[part->indices[i]] = current_side;
            std::cout << "\ttrue " << part->indices[i] << std::endl;
        }
    // Sort vertices (x axis)
        for (int i = 0; i < part->indices.size(); i ++)
        {
            events.push_back(SubPolygon::Vertex(i, &(*part)));
        }
        std::sort(events.begin(), events.end(), less_than<X, SubPolygon::Vertex>);

    //TRIANGULATE the subpolygon 'part'
        std::list<SubPolygon::Vertex> L;
        L.push_back(events[0]); L.push_back(events[1]);
        for (auto vertex = events.begin() + 2; vertex != events.end(); vertex ++)
        {
            std::cout << "-* " << vertex->getIndex() << "*-" << std::endl;

            if (side[vertex->getIndex()] == side[L.back().getIndex()])
            {
                // As long as angle between vertex and two last vertices in list is convex..
                SubPolygon::Vertex second_last = *(++L.rbegin());
                while (L.size() > 1
                        && leftof( *L.back() - *second_last, **vertex - *L.back()) ^ side[vertex->getIndex()])
                {
                    second_last = *(++ L.rbegin());
                    // Add triangle (last, next-to-last, vertex)
                    // diagonals.push_back( Diagonal(vertex->getIndex(), second_last.getIndex(), this));
                    float f = randFloat();
                    triangles.push_back( Triangle(**vertex, *second_last, *L.back(), f * color));
                    // std::cout << "Triangle(" << vertex->getIndex() << ", " << second_last.getIndex() << ", " << L.back().getIndex() << ")" << std::endl;
                    L.pop_back();
                }
                L.push_back(*vertex);
            }
            else // The vertex is on the other side of vertices in L
            {
                SubPolygon::Vertex second = *(++ L.begin());
                while (L.size() > 1)
                {
                    second = *(++L.begin());
                    // Add triangle (first, second, vertex)
                    // diagonals.push_back( Diagonal(vertex->getIndex(), second.getIndex(), this));
                    triangles.push_back( Triangle(**vertex, *second, *L.front(), randFloat() * color));
                    // std::cout << "Triangle(" << vertex->getIndex() << ", " << second.getIndex() << ", " << L.front().getIndex() << ")" << std::endl;
                    L.pop_front();
                }
                L.push_back(*vertex);
            }
        }
    }
}

int Polygon::numEdges()
{
	return vertices.size();
}
//TODO Caution: Local coordinates
LineSegment Polygon::getEdge(int index)
{
	if (index < 0 || index >= numEdges()) throw(std::out_of_range("Edge out of range."));
	glm::vec2 s = vertices[index];
	glm::vec2 t = vertices[(index < vertices.size() - 1)? index + 1 : 0];
	// s += position;
	// t += position;
	return LineSegment(s, t);
}

//////////////////////////////////////////////
//
//              Vertex pointer
//
//////////////////////////////////////////////

Polygon::Vertex::Vertex(int index, Polygon* parent)
    :index(index % parent->vertices.size()), parent(parent) { }

glm::vec2& Polygon::Vertex::operator* ()
{
	return parent->vertices[index];
}
glm::vec2* Polygon::Vertex::operator-> ()
{
	return &(parent->vertices[index]);
}
glm::vec2& Polygon::Vertex::preceding()
{
	int i = (index == 0) ? parent->vertices.size() - 1 : index - 1;
	return *(&(parent->vertices[i]));
}
glm::vec2& Polygon::Vertex::successive()
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}

//////////////////////////////////////////////
//
//              Edge pointer
//
//////////////////////////////////////////////

Polygon::Edge::Edge(int index, Polygon *parent)
    :index(index % parent->vertices.size()), parent(parent) { }

glm::vec2& Polygon::Edge::start() const
{
	return parent->vertices[index];
}
glm::vec2& Polygon::Edge::end() const
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}

// Find y value given x value. Doesn't care about bounds.
int Polygon::Edge::operator() (int x) const
{
	glm::vec2 delta = end() - start();
	assert(start().x != end().x);
	float t = (x - start().x)/delta.x;
	return start().y + t * delta.y;
}

bool Polygon::Edge::operator== (Edge other)
{
    return index == other.index && parent == other.parent;
}


//////////////////////////////////////////////
//
//              SubPolygon
//
//////////////////////////////////////////////
//

void SubPolygon::fill()
{
    indices.clear();
    for (int i = 0; i < mother->vertices.size(); i ++) {
        indices.push_back(i);
    }
}
bool SubPolygon::containsDiagonal(int a, int b)
{
    int size = indices.size();
    int a_index = -1, b_index = -1;
    for (int i = 0; i < size; i ++)
    {
        if (indices[i] == a) a_index = i;
        if (indices[i] == b) b_index = i;
    }
    bool neighbors = ((a_index + 1) % size == b_index) || ((b_index + 1) % size == a_index);
    if (a_index != -1 && b_index != -1 && !neighbors) {
        return true;
    }
    return false;
}

void SubPolygon::split(int a, int b, SubPolygon& out1, SubPolygon& out2)
{
    out1.indices.clear();
    out2.indices.clear();
    // Find indices
    int i_a = -1, i_b = -1;
    for (int i = 0; i < indices.size(); i ++) {
        if (indices[i] == a) i_a = i;
        if (indices[i] == b) i_b = i;
    }
    assert(i_a != -1 && i_b != -1);
    // range [i_a, i_b]
    for (int i = i_a; ; i ++) {
        if (i >= indices.size()) i = 0;
        out1.indices.push_back(indices[i]);
        if (i == i_b) break;
    }
    // range [i_b, i_a]
    for (int i = i_b; ; i ++) {
        if (i >= indices.size()) i = 0;
        out2.indices.push_back(indices[i]);
        if (i == i_a) break;
    }
}

std::ostream& operator<< (std::ostream& lhs, SubPolygon& p)
{
    lhs << "SubPolygon(";
    for (auto it = p.indices.begin(); it != p.indices.end(); it ++)
    {
        lhs << *it << ", ";
    }
    lhs << ")\n";
    return lhs;
}


///////////////////////
// Vertex smart-pointer
///////////////////////

SubPolygon::Vertex::Vertex(int index, SubPolygon *parent)
    :index(index), parent(parent) { }


glm::vec2& SubPolygon::Vertex::operator* ()
{
	return parent->mother->vertices[parent->indices[index]];
}
glm::vec2* SubPolygon::Vertex::operator-> ()
{
	return &(parent->mother->vertices[parent->indices[index]]);
}
glm::vec2& SubPolygon::Vertex::preceding()
{
	int i = (index == 0) ? parent->indices.size() - 1 : index - 1;
	return *(&(parent->mother->vertices[parent->indices[i]]));
}
glm::vec2& SubPolygon::Vertex::successive()
{
	int i = (index == parent->indices.size() - 1) ? 0 : index + 1;
	return parent->mother->vertices[parent->indices[i]];
}

///////////////////////
// Edge smart-pointer
///////////////////////

SubPolygon::Edge::Edge(int index, SubPolygon *parent)
    :index(index % parent->indices.size()), parent(parent) { 
}

glm::vec2& SubPolygon::Edge::start() const
{
	return parent->mother->vertices[parent->indices[index]];
}
glm::vec2& SubPolygon::Edge::end() const
{
	int i = (index == parent->indices.size() - 1) ? 0 : index + 1;
	return parent->mother->vertices[parent->indices[i]];
}

int SubPolygon::Edge::getIndex() const
{
    return parent->indices[index]; 
}

// Find y value given x value. Doesn't care about bounds.
int SubPolygon::Edge::y(int x) const
{
	glm::vec2 delta = end() - start();
    if (delta.x == 0) {
        return delta.y / 2;
    }
	float t = (x - start().x)/delta.x;
	return start().y + t * delta.y;
}
int SubPolygon::Edge::x(int y) const
{
	glm::vec2 delta = end() - start();
    if (delta.y == 0) {
        return delta.x / 2;
    }
	assert(start().y != end().y);
	float t = (y - start().y)/delta.y;
	return start().x + t * delta.x;
}

bool SubPolygon::Edge::operator== (Edge other) const 
{
    return index == other.index && parent == other.parent;
}
