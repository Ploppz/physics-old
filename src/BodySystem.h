#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <iostream>

#include "Polygon.h"
#include "glutils.h"

using namespace glm;
using namespace std;

class Body;

class BodySystem
{
	friend class Body;
public:
	BodySystem();

	void timestep(float delta);
	Body addBody(); // Returns index of new body
    int numBodies() const { return count; }
	Body getBody(int index);
	Polygon& getPolygon(int body);

	bool overlaps(Body a, Body b);
private:
	int count;

	vector<float> mass;

	vector<vec2> position;
	vector<vec2> velocity;
	vector<vec2> force;

	vector<float> orientation;
	vector<float> angularSpeed;
	vector<float> torque;

	// Temporary solution: only one polygon each Body
	vector<Polygon> shape;
};

// class BODY
// Functions as an index
// Flaws: Gets invalid if elements move within vector
class Body
{
	friend class BodySystem;
public:

	// Should invalidate it.
	Body() {}

	// Access to members
		float& mass() { return system->mass[index]; }
		vec2& position() { return system->position[index]; }
		vec2& velocity() { return system->velocity[index]; }
		vec2& force() { return system->force[index]; }

		float& orientation() { return system->orientation[index]; }
		float& angularSpeed() { return system->angularSpeed[index]; }
		float& torque() { return system->torque[index]; }

		Polygon& shape() { return system->shape[index]; }
	// Drawing
        void addToBuffer(BufferWriter<float> &buffer);
		void addToBuffer(std::vector<float> &buffer);
		unsigned int addToBuffer(float *buffer, int offset);
	// Iterator TODO

private:
	Body(BodySystem *system, int index);

	BodySystem *system;
	int index;
};
