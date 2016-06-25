#pragma once
#include "../geometry/EdgePoint.h"
struct DepthContact {
    float depth;
    glm::vec2 normal;
    EdgePoint subj_point;
    EdgePoint ref_point;
};
// rewind_time is positive for 'back in time'
struct TimeContact {
	float rewind_time;
	// Collision normal in world coordinates
	glm::vec2 normal;

	// Collision incident points along the polygons
	EdgePoint ref_point;
	EdgePoint subj_point;
	TimeContact() {};
	TimeContact(DepthContact other) : rewind_time(0), normal(other.normal), ref_point(other.ref_point), subj_point(other.subj_point) {};
};
