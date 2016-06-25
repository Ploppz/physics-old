#pragma once
// Test arena for physics
// Test case:
// no friction, no gravity, bounded by window size

#include <vector>
#include "BodySystem.h"

class World 
{
public:

	BodySystem bodies;

	void timestep(float delta)
	{
		bodies.timestep(delta);
		// TODO: Make collision test & resolution for boundaries

	}
};
