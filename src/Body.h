#pragma once
#include "geometry/Polygon.h"
#include "geometry/EdgePoint.h"
#include "BodySystem.h"


/** Body -- Smartpointer **/

// These should maybe be defined in BodySystem.h...

class Body
{
	friend class BodySystem;
 public:
		PolygonMode mode;

	// Should invalidate it.
		Body() {mode = POLYGON_INSIDE;}
	// General
		void apply_impulse(glm::vec2 impulse, EdgePoint point);
		void apply_impulse(glm::vec2 impulse, glm::vec2 point); //TODO implement? don't know if I need it

        void update(float delta_time);
		void update_polygon_state();

	// Access to members
		glm::vec2 real_position(); 
		/* inline float& mass() { return system->mass[index]; } */
		inline PositionType& position_type() { ensure_valid(); return system->position_type[index]; }
		inline glm::vec2& position() { ensure_valid(); return system->position[index]; }
		inline glm::vec2& velocity() { ensure_valid(); return system->velocity[index]; }
		inline glm::vec2& force() { ensure_valid(); return system->force[index]; }

		inline float& orientation() { ensure_valid(); return system->orientation[index]; }
		inline float& rotation() { ensure_valid(); return system->rotation[index]; }
		inline float& torque() { ensure_valid(); return system->torque[index]; }

		inline Polygon& shape() { ensure_valid(); return system->shape[index]; }

		inline Body& parent() { ensure_valid(); return system->parent[index]; }

	// Iterator?
		Body& operator++ ();
		bool is_valid();
		int get_index() {return index; }

 private:
	Body(int index, BodySystem *system);
	void ensure_valid();

	BodySystem *system;
	int index;
};
std::ostream& operator<< (std::ostream& out, Body b);
