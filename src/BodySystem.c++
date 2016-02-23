
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <glm/gtx/matrix_transform_2d.hpp>

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
        shape[i].matrix = glm::translate(glm::mat3 {}, position[i]);
	}
}

Body BodySystem::addBody()
{
    position_type.push_back(RELATIVE);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
	angularSpeed.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon () );

    children.push_back(std::vector<Body> {});
    parent.push_back(Body(this, -1)); // Invalid body

	++ count;
	return Body(this, count - 1);
}
Body BodySystem::addBody(Body parent)
{

    position_type.push_back(RELATIVE);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
	angularSpeed.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon () );

    assert(parent.index <= int(children.size() - 1));
    children[parent.index].push_back(Body(this, count));
    children.push_back(std::vector<Body> {});
    this->parent.push_back(parent);

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


// class BODY

Body::Body(BodySystem *system, int index)
	: system(system), index(index)
{
	// std::cout << "Body .. " << index << std::endl;
}

vec2 Body::real_position()
{
    if (position_type() == ABSOLUTE || parent().index == -1) {
        return position();
    } else { // RELATIVE
        return position() + parent().real_position();
    }
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
		// buffer.push_back(shape().color.r);
		// buffer.push_back(shape().color.g);
		// buffer.push_back(shape().color.b);
		// Draw to next point
		if (next != shape().vertices.end()) {
			buffer.push_back(next->x);
			buffer.push_back(next->y);
		} else {
			// End the shape.
			buffer.push_back(first->x);
			buffer.push_back(first->y);
		}
		// buffer.push_back(shape().color.r);
		// buffer.push_back(shape().color.g);
		// buffer.push_back(shape().color.b);
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
		// buffer.write(shape().color.r, shape().color.g, shape().color.b);
		// Draw to next point
		if (next != shape().vertices.end()) {
            buffer.write(next->x, next->y);
		} else {
			// End the shape.
			buffer.write(first->x, first->y);
		}
		// buffer.write(shape().color.r, shape().color.g, shape().color.b);

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
