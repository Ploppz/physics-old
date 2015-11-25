
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>

#include "BodySystem.h"
#include "glutils.h"
#include "typewriter/FontRenderer.h"
#include "Geometry.h"

extern FontRenderer *fontRenderer;

using namespace glm;

BodySystem::BodySystem()
	: count {}, mass {}, position {}, velocity {}, force {}, orientation {}, angularSpeed {}, torque {}, shape {}
{
}

void BodySystem::timestep(float delta)
{
	for (int i = 0; i < count; i ++)
	{
		position[i] += velocity[i];
	}
}

Body BodySystem::addBody()
{
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
	angularSpeed.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon () );

	++ count;
	return Body(this, count - 1);
}

Body BodySystem::getBody(int index)
{
	return Body(this, index);
}

// Collision

typedef std::pair<glm::vec2, glm::vec2> LineSegment;

bool inside(vec2 point, Polygon polygon);

// For now this is copy-pasta code from Polygon: we just have to add the Body's pos.
bool BodySystem::overlaps(Body a, Body b)
{
	// Make an "Influence Area" from a
	// test centroid of b.
	glm::vec2 centroid_a = a.shape().centroid() + a.position();
	glm::vec2 centroid_b = b.shape().centroid() + b.position();
	float radius_b = b.shape().radius();
	// First, loop through edges (u, v) of a, and create a triangle with the centroid, from which we
	// calculate the barycentric coordinate space
    

    if (inside(b.shape().vertices[0] + b.position(), a)) {
        return true;
    }

	float distFromEdge;
	float slope, minVal, deltaVal;
	glm::vec2 delta;
	bool withinBounds;
	LineSegment edge_a, edge_b;
	for (int i = 0; i < a.shape().numEdges(); i ++)
	{
		edge_a = a.shape().getEdge(i);
		edge_a.first += a.position(); edge_a.second += a.position();
		// Get barycentric coordinates with respect to centroid_a
		glm::vec3 lala = barycentric(centroid_a, edge_a.first, edge_a.second, centroid_b);
		distFromEdge = distance(centroid_b, edge_a.first, edge_a.second);

		// If absolute value of this is within the radius of b, then there is a possible collision
		
		{	// First, also limit the y or x coordinate based on slope
			delta = edge_a.first - edge_a.second;
			slope = delta.y / delta.x;
			if (abs(slope) < 1) {	// Limit on x axis
				minVal = std::min(edge_a.first.x, edge_a.second.x) - radius_b;
				deltaVal = abs(delta.x) + radius_b + radius_b;
				withinBounds = centroid_b.x > minVal && centroid_b.x < minVal + deltaVal;
			} else {				// Limit on y axis
				minVal = std::min(edge_a.first.y, edge_a.second.y) - radius_b;
				deltaVal = abs(delta.y) + radius_b + radius_b;
				withinBounds = centroid_b.y > minVal && centroid_b.y < minVal + deltaVal;
			}
		}
		if (abs(distFromEdge) <= radius_b && withinBounds) {
			// DO DETAILED TEST b vs. this edge (s, t)
			//TODO: eliminate edges where both vertices are outside (using the bary-space of the centroid-triangle)
			//TODO: ... we already know which edges may collide and which not
			// For now: exhaustive. Loop through all edges of b.
			for (int j = 0; j < b.shape().numEdges(); j ++)
			{
				edge_b = b.shape().getEdge(j);
				edge_b.first += b.position(); edge_b.second += b.position();
				if (intersect(edge_a.first, edge_a.second, edge_b.first, edge_b.second)) {
					return true;
				}
			}
		}
	}
	return false;
}






// class BODY

Body::Body(BodySystem *system, int index)
	: system(system), index(index)
{
	// std::cout << "Body .. " << index << std::endl;
}






// DRAWING

using namespace std;

void Body::addToBuffer(std::vector<float> &buffer)
{
	// Test: render lines.
	std::vector<glm::vec2>::iterator next;
    int i = 0;
	auto first = shape().vertices.begin();
	for (auto it = first; it != shape().vertices.end(); it ++)
    {
		next = it; next ++;
		buffer.push_back(it->x);
		buffer.push_back(it->y);
		buffer.push_back(shape().color.r);
		buffer.push_back(shape().color.g);
		buffer.push_back(shape().color.b);
		// Draw to next point
		if (next != shape().vertices.end()) {
			buffer.push_back(next->x);
			buffer.push_back(next->y);
		} else {
			// End the shape.
			buffer.push_back(first->x);
			buffer.push_back(first->y);
		}
		buffer.push_back(shape().color.r);
		buffer.push_back(shape().color.g);
		buffer.push_back(shape().color.b);
        // Add text
        fontRenderer->addText(static_cast<ostringstream*>( &(ostringstream() << i) )->str(), position().x + it->x, position().y + it->y,  false);
        i ++;
	}
	
	// RENDER RADIUS
	/* float r = shape().radius();
	glm::vec2 c = shape().centroid() + position();
	for (int i = 0; i < 50; i ++) {
		buffer.push_back(c.x + cos(i / 25.0f * M_PI) * r);
		buffer.push_back(c.y + sin(i / 25.0f * M_PI) * r);
		// Color
		buffer.push_back(0);
		buffer.push_back(1);
		buffer.push_back(0);
		//
		buffer.push_back(c.x + cos((i + 1) / 25.0f * M_PI) * r);
		buffer.push_back(c.y + sin((i + 1) / 25.0f * M_PI) * r);
		// Color
		buffer.push_back(0);
		buffer.push_back(1);
		buffer.push_back(0);
	} */
}
void Body::addToBuffer(BufferWriter<float> &buffer)
{
	// Test: render lines.
	std::vector<glm::vec2>::iterator next;
    int i = 0;
	auto first = shape().vertices.begin();
	for (auto it = first; it != shape().vertices.end(); it ++)
    {
		next = it; next ++;
		buffer.write(it->x, it->y);
		buffer.write(shape().color.r, shape().color.g, shape().color.b);
		// Draw to next point
		if (next != shape().vertices.end()) {
            buffer.write(next->x, next->y);
		} else {
			// End the shape.
			buffer.write(first->x, first->y);
		}
		buffer.write(shape().color.r, shape().color.g, shape().color.b);

        // Add text
        fontRenderer->addText(static_cast<ostringstream*>( &(ostringstream() << i) )->str(), position().x + it->x, position().y + it->y,  false);
        i ++;
	}
}
unsigned int Body::addToBuffer(float *buffer, int offset)
{
	std::vector<float> data;
	addToBuffer(data);
	// Write data to buffer
	
	std::copy(data.begin(), data.end(), buffer + offset);
	return offset + data.size();
}
