#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <iostream>

#include "geometry/Polygon.h"
#include "glutils.h"
#include "geometry/geometry.h"


enum PositionType { ABSOLUTE, RELATIVE };

class Body;

class BodySystem
{
	friend class Body;
public:
	BodySystem();

	void timestep(float delta);
	Body add_body(); // Returns index of new body
    Body add_body(Body parent);
    int num_bodies() const { return count; }
	Body get_body(int index);
	Polygon& get_polygon(int body);

    std::vector<Intersect> overlaps(Body a, Body b);
private:
	int count;

    // Body data
	std::vector<float> mass;

    std::vector<PositionType> position_type;

	std::vector<glm::vec2> position;
	std::vector<glm::vec2> velocity;
	std::vector<glm::vec2> force;

	std::vector<float> orientation;
	std::vector<float> angular_speed;
	std::vector<float> torque;

	// Temporary solution: only one polygon each Body
	std::vector<Polygon> shape;

    // Hierarchical tree of bodies
    std::vector<Body> parent;
    std::vector<std::vector<Body>> children;
};

// class BODY
// Functions as an index
// Flaws: Gets invalid if elements move within std::vector
class Body
{
	friend class BodySystem;
public:

	// Should invalidate it.
	Body() {}

	// Access to members
		inline float& mass() { return system->mass[index]; }
        inline PositionType& position_type() { return system->position_type[index]; }
		inline glm::vec2& position() { return system->position[index]; }
        glm::vec2 real_position();
		inline glm::vec2& velocity() { return system->velocity[index]; }
		inline glm::vec2& force() { return system->force[index]; }

		inline float& orientation() { return system->orientation[index]; }
		inline float& angular_speed() { return system->angular_speed[index]; }
		inline float& torque() { return system->torque[index]; }

		inline Polygon& shape() { return system->shape[index]; }

        inline Body& parent() { return system->parent[index]; }

	// Drawing
        void add_to_buffer(BufferWriter<float> &buffer);
		void add_to_buffer(std::vector<float> &buffer);
		unsigned int add_to_buffer(float *buffer, int offset);
	// Iterator TODO

private:
	Body(BodySystem *system, int index);

	BodySystem *system;
	int index;
};
