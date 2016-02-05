#include "Polygon.h"
#include "Geometry.h"
#include "tmp.h"
#include "LinkedList.h"
#include "glutils.h"

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
    glm::vec2 a, b;
	std::vector<glm::vec2>::iterator next;
	for (auto it = vertices.begin(); it != vertices.end(); it ++)
	{
		next = it; next ++;
		if (next == vertices.end()) break;
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

glm::vec2 Polygon::transform(glm::vec2 point)
{
    return glm::vec2(   matrix[0][0] * point.x + matrix[1][0] * point.y + matrix[2][0],
                        matrix[0][1] * point.x + matrix[1][1] * point.y + matrix[2][1] );
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
	if (rhs.start().x < rhs.end().x)    b_leftmost = rhs.start();
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
            triangulate(parts, diagonals, triangles);
        }

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

	for (int i = 0; i < events.size() - 1; i ++) // Each 'event[i]' is a vertex.
	{

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
		
        if (precSide == succSide  &&  dir != precSide) // Both adjacent vertices are BEHIND the scanline
		{
			if (concave) { /* MERGE POINT */
				// Remove the top-most edge from status
				index = events[i].getIndex() - (logic);
				if (index == -1) index = indices.size() - 1;
                erase(status, index, this);
				// Update helper of the edge that is just below this vertex
                Edge justBelow = getEdgeJustBelow<axis>(status, events[i]);
				helpers[justBelow.getIndex()] = events[i];
			} else { /* END POINT */
				// Remove the bottom-most edge from status
				index = events[i].getIndex() - (logic);
				if (index == -1) index = indices.size() - 1;
                erase(status, index, this);
			}
		} else if (precSide == succSide  &&  dir == precSide) // Both adjacent vertices are IN FRON OF the scanline
		{
			if (concave) { /* SPLIT POINT */
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
				// Insert the bottom-most incident edge to 'status'
				index = events[i].getIndex() - !(logic);
				if (index == -1) index = indices.size() - 1;
                Edge e = Edge(index, this);
				status.insert( e);
				helpers[index] = events[i]; // with this vertex as helper
			}
		} else { /* NORMAL VERTEX */
			// CW: If the incident edge (any of the two) goes to the right, then this vertex is in the _ceiling_
			bool ceiling;
            ceiling = (precSide == BACK && succSide == FORTH) ^ CCW;
            if (axis == X) {
                // ceiling = !(events[i]->x > events[i].successive().x) ^ CCW;
            } else {
                // ceiling = !(events[i]->y < events[i].successive().y) ^ CCW;
            }

			if (ceiling) {
				// Replace helper of the edge just below
                Edge justBelow = getEdgeJustBelow<axis>(status, events[i]);
				helpers[justBelow.getIndex()] = events[i];
			} else {
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
        for (uint i = min_index; ; i ++) {
            if (i >= part->indices.size()) i = 0;
            if (i == max_index) break;
            side[part->indices[i]] = current_side;
        }
        current_side = !current_side;
        for (int i = max_index; ; i ++) {
            if (i >= part->indices.size()) i = 0;
            if (i == min_index) break;
            side[part->indices[i]] = current_side;
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
                    float intensity = 0.5 + randFloat()*0.5;
                    triangles.push_back( Triangle(**vertex, *second_last, *L.back(), glm::vec3 {}));
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
                    float intensity = 0.5 + randFloat()*0.5;
                    triangles.push_back( Triangle(**vertex, *second, *L.front(), glm::vec3 {}));
                    // std::cout << "Triangle(" << vertex->getIndex() << ", " << second.getIndex() << ", " << L.front().getIndex() << ")" << std::endl;
                    L.pop_front();
                }
                L.push_back(*vertex);
            }
        }
    }
}

std::vector<Intersection> Polygon::overlaps(Polygon& a, Polygon& b)
{
    bool collision = false;
	// Make an "Influence Area" from a
	// test centroid of b.
	glm::vec2 centroid_a = a.transform(a.centroid());
	glm::vec2 centroid_b = b.transform(b.centroid());
	float radius_b = b.radius();
	// First, loop through edges (u, v) of a, and create a triangle with the centroid, from which we
	// calculate the barycentric coordinate space
    

    if (inside(b.transform(b.vertices[0]), a) ||
        inside(a.transform(a.vertices[0]), b)) {
        collision = true;
    }
    std::vector<Intersection> intersections;

	float distFromEdge;
	float slope, min_val, delta_val;
	glm::vec2 delta;
	bool withinBounds;
	LineSegment edge_a, edge_b;
	for (int i = 0; i < a.numEdges(); i ++)
	{
		edge_a = a.getEdge(i);
        edge_a.first = a.transform(edge_a.first);
        edge_a.second = a.transform(edge_a.second);
		// Get barycentric coordinates with respect to centroid_a
		// glm::vec3 lala = barycentric(centroid_a, edge_a.first, edge_a.second, centroid_b);

		// If absolute value of this is within the radius of b, then there is a possible collision
		
		{	// First, also limit the y or x coordinate based on slope
			delta = edge_a.first - edge_a.second;
			slope = delta.y / delta.x;
			if (fabs(slope) < 1) {	// Limit on x axis
				min_val = std::min(edge_a.first.x, edge_a.second.x) - radius_b;
				delta_val = fabs(delta.x) + radius_b + radius_b;
				withinBounds = centroid_b.x > min_val && centroid_b.x < min_val + delta_val;
			} else {				// Limit on y axis
				min_val = std::min(edge_a.first.y, edge_a.second.y) - radius_b;
				delta_val = fabs(delta.y) + radius_b + radius_b;
				withinBounds = centroid_b.y > min_val && centroid_b.y < min_val + delta_val;
			}
		}
		distFromEdge = distance(centroid_b, edge_a.first, edge_a.second);
		if (fabs(distFromEdge) <= radius_b && withinBounds) {
			// DO DETAILED TEST b vs. this edge (s, t)
			//TODO: eliminate edges where both vertices are outside (using the bary-space of the centroid-triangle)
			//TODO: ... we already know which edges may collide and which not
			for (int j = 0; j < b.numEdges(); j ++)
			{
				edge_b = b.getEdge(j);
                edge_b.first = b.transform(edge_b.first);
                edge_b.second = b.transform(edge_b.second);
				if (intersect(edge_a.first, edge_a.second, edge_b.first, edge_b.second)) {
                    intersections.push_back(Intersection(   Polygon::Edge(i, &a),
                                                            Polygon::Edge(j, &b)));
                    collision = true;
				}
			}
		}
	}
	return intersections;
}

// INTERSECTION

typedef bool Side;
const bool IN = true;
const bool OUT = false;
struct NewVertex {
    NewVertex(Intersection& i, Side in_out) {
        coor = i.point;
        intersect = true;
        processed = false;
        this->in_out = in_out;
    }
    NewVertex(Polygon::Vertex vert, Side in_out) {
        coor = *vert;
        intersect = false;
        this->in_out = in_out;
    }

    glm::vec2 coor;
    bool intersect; // Intersection point or normal vertex
    // If intersect:
        bool processed;
    Side in_out; // Entry or exit to other polygon?

    NewVertex *prev, *next, *parallel; // Linked list, and link to evt. same intersection point in other polygon
};

struct FullIntersect { // An intersection is needed to sort the NewVertices
    FullIntersect(Intersection *i, NewVertex *p, NewVertex *q): i(i), vert_p(p), vert_q(q){}
    Intersection* i;
    NewVertex* vert_p;
    NewVertex* vert_q;
    template <Intersection::Which which> // sort wrt p or q?
    static bool lt(FullIntersect& a, FullIntersect& b);
};


////////////////////////////
// CALCULATE INTERSECTION //
///////////////////////////

std::vector<Polygon> Polygon::intersection(Polygon& p, Polygon& q)
{
    std::vector<Intersection> intersects = overlaps(p, q); // Let an intersect be an intersection vertex

    /** LINK PARALLEL VERTICES **/
    std::vector<FullIntersect> sorted;
    for (auto it = intersects.begin(); it != intersects.end(); it ++)
    {
        Intersection i = *it;
        NewVertex *p = new NewVertex(i, OUT);
        NewVertex *q = new NewVertex(i, OUT);
        p->parallel = q; q->parallel = p;
        sorted.push_back(  FullIntersect(&i, p, q)  );
    }

    /** Sort wrt p **/
    //  The FIRST polygon of Intersection is from Polygon p
    std::sort(sorted.begin(), sorted.end(), FullIntersect::lt<Intersection::FIRST>);
    /** Interleave vertices and intersects in p **/
    CircularList<NewVertex> p_vertices;

    Side in_out = inside(p.vertices[0], q);
    int added_vertex_index = -1;
    int current_index = 0;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge1.getIndex();
        // Add eventual vertices that were skipped
        Vertex vertex(current_index, &p);
        while (current_index > added_vertex_index) {
            p_vertices.push_back(new NewVertex(vertex, in_out));

            ++ vertex;
            ++ added_vertex_index;
        }
        // Add the intersection point
        in_out = !in_out;
        it->vert_p->in_out = in_out;
        p_vertices.push_back(it->vert_p);
    }

    /** Sort wrt q **/
    //  The FIRST polygon of Intersection is from Polygon p
    std::sort(sorted.begin(), sorted.end(), FullIntersect::lt<Intersection::SECOND>);
    /** Interleave vertices and intersects in q **/
    CircularList<NewVertex> q_vertices;

    in_out = inside(q.vertices[0], q);
    added_vertex_index = -1;
    current_index = 0;
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        current_index = it->i->edge1.getIndex();
        // Add eventual vertices that were skipped
        Vertex vertex(current_index, &q);
        while (current_index > added_vertex_index) {
            q_vertices.push_back(new NewVertex(vertex, in_out));

            ++ vertex;
            ++ added_vertex_index;
        }
        // Add the intersection point
        in_out = !in_out;
        it->vert_q->in_out = in_out;
        q_vertices.push_back(it->vert_q);
    }

    std::vector<Polygon> result;
    /** Create Polygons - by traversing q (because we already have the intersection points of q sorted. **/
    for (auto it = sorted.begin(); it != sorted.end(); it ++)
    {
        Polygon polygon;
       
        NewVertex* start;
        NewVertex* current = start = it->vert_q;
        // current and start are _intersection_
        //std::cout << std::boolalpha << current->processed << std::endl;
        
        if (current->processed) continue;

        //std::cout << "LA" << std::endl;
        Direction direction = current->in_out;

        // This is the start of a new polygon.
        do {
            polygon.vertices.push_back(current->coor);

            if (current->intersect) {
                assert(current->parallel);
                current->processed = true;
                current->parallel->processed = true;
                current = current->parallel;
                direction = current->in_out;
            }

            if (direction == FORTH) {
                assert(current->next);
                current = current->next;
            } else {
                assert(current->prev);
                current = current->prev;
            }
        } while (current != start);


        result.push_back(polygon);
    }

    return result;
}


// TRIANGLE FAN TODO make this clear
//void Polygon::appendStencilTriangles(std::vector<float> &buffer)
void Polygon::appendStencilTriangles(BufferWriter<float> &buffer)
{
    // i is the edge index
    for (int i = 0; i < vertices.size(); i ++)
    {
        buffer.write(vertices[i].x, vertices[i].y);
    }
}
void Polygon::appendLinesToVector(std::vector<float> &list)
{
    for (int i = 0; i < vertices.size(); i ++) {
        int j = i + 1; j %= vertices.size();
        glm::vec2 vec_i = transform(vertices[i]);
        glm::vec2 vec_j = transform(vertices[j]);
        list.push_back(vec_i.x);
        list.push_back(vec_i.y);
        list.push_back(vec_j.x);
        list.push_back(vec_j.y);
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
glm::vec2& Polygon::Vertex::successive()
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->vertices[i];
}
glm::vec2& Polygon::Vertex::preceding()
{
	int i = (index == 0) ? parent->vertices.size() - 1 : index - 1;
	return *(&(parent->vertices[i]));
}
Polygon::Vertex& Polygon::Vertex::operator++ ()
{
	index = (index == parent->vertices.size() - 1) ? 0 : index + 1;
    return *this;
}
Polygon::Vertex& Polygon::Vertex::operator-- ()
{
	index = (index == 0) ? parent->vertices.size() - 1 : index - 1;
    return *this;
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
glm::vec2 Polygon::Edge::start_tr() const
{
	return parent->transform(parent->vertices[index]);
}
glm::vec2 Polygon::Edge::end_tr() const
{
	int i = (index == parent->vertices.size() - 1) ? 0 : index + 1;
	return parent->transform(parent->vertices[i]);
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
// Vertex (SubPolygon)
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
// Edge (SubPolygon)
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




////////////////////
// Intersection
///////////////////
Intersection::Intersection(Polygon::Edge edge1, Polygon::Edge edge2)
{
    this->edge1 = edge1; this->edge2 = edge2;
    // Find intersection point (last argument is output)
    bool result = intersect(edge1.start_tr(), edge1.end_tr(),
            edge2.start_tr(), edge2.end_tr(), point, alpha1, alpha2);
    if (!result) {
        std::cout << "Error: edges do not actually overlap (" << edge1.getIndex() << ", " << edge2.getIndex() << ")" << std::endl;
    }
}

template <Intersection::Which which> // TODO separate functions for FIRST and SECOND
bool Intersection::lt(Intersection& i, Intersection& j)
{
    if (which == FIRST) {
        if (i.edge1.getIndex() < j.edge1.getIndex()) {
            return true;
        } else if (i.edge1.getIndex() == j.edge1.getIndex()) {
            return (i.alpha1 < j.alpha1);
        } else {
            return false;
        }
    } else {
        if (i.edge2.getIndex() < j.edge2.getIndex()) {
            return true;
        } else if (i.edge2.getIndex() == j.edge2.getIndex()) {
            return (i.alpha2 < j.alpha2);
        } else {
            return false;
        }
    }
}

template <Intersection::Which which>
bool FullIntersect::lt(FullIntersect& a, FullIntersect& b)
{
    Intersection::lt<which>(*(a.i), *(b.i));
}
