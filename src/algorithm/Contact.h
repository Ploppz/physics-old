#pragma once
#include "geometry/EdgePoint.h"

// Notes:
// It's kind of unnecessary to use inheritance here.
// E.g. TimeContact doesn't need make_normal_point_toward_subject
//
struct DepthContact {
    float depth;
    glm::vec2 normal;
    EdgePoint subj_point;
    EdgePoint ref_point;

    void ensure_normal_toward_subject();

};
// rewind_time is positive for 'back in time'
struct TimeContact {
	TimeContact() {};
	TimeContact(DepthContact other)
        : rewind_time(0), normal(other.normal), ref_point(other.ref_point), subj_point(other.subj_point) {};

	float rewind_time;
	glm::vec2 normal;
	EdgePoint ref_point;
	EdgePoint subj_point;
};


/*
Disregards time, finds the distance between the edge & vertex incident.
*/
DepthContact simple_contact_conversion(TimeContact& in_contact);
