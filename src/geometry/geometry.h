#pragma once

#include <glm/glm.hpp>
#include <iostream>

#include "../BodySystem.h"
#include "Polygon.h"

class Body;

std::ostream &operator << (std::ostream &lhs, const glm::vec2 &rhs);
std::ostream &operator << (std::ostream &lhs, const glm::vec3 &rhs);

float angle(glm::vec2 a, glm::vec2 b, glm::vec2 c);

inline glm::vec2 unit_vector(float angle) {
    return glm::vec2(cos(angle), sin(angle));
}
inline glm::vec2 unit_vector_wrt_y(float angle) {
    return glm::vec2(- sin(angle), cos(angle));
}
inline float angle_of_vector(glm::vec2 vec) {
    return atan2(vec.y, vec.x);
}
// Returns true if a points left of b
bool leftof(glm::vec2 a, glm::vec2 b);

float distance(glm::vec2 a, glm::vec2 b);
float distance_line(glm::vec2 p, glm::vec2 line_a, glm::vec2 line_b);

// Distance from polygon (also gives the closest edge)
float distance(glm::vec2 point, Polygon& p, int& out_closest_edge, float& out_closest_edge_alpha);
float distance(glm::vec2 point, Polygon& p);
float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end, float &out_alpha);
float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end);

float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end, float &out_alpha);

// Midpoint of two points
glm::vec2 middle(glm::vec2 a, glm::vec2 b);

/* Projection */
glm::vec2 project(glm::vec2 p, glm::vec2 line_a, glm::vec2 line_b);
float project(glm::vec2 point, glm::vec2 direction);


float signedArea(glm::vec2 a, glm::vec2 b, glm::vec2 c);

// First 3 points are for triangle.
glm::vec3 barycentric(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 p);
bool intersect(glm::vec2 line1_a, glm::vec2 line1_b, glm::vec2 line2_a, glm::vec2 line2_b);
bool intersect(glm::vec2 line1_a, glm::vec2 line1_b, glm::vec2 line2_a, glm::vec2 line2_b, glm::vec2& point_of_intersection, float& alpha1, float& alpha2);
// return the x value of the intersection between a line and a constant y value
float intersect_horizontal(glm::vec2 line_start, glm::vec2 line_direction, float y_constant, float &alpha_out);

bool inside(glm::vec2 point, Polygon& p);
float cross(glm::vec2 a, glm::vec2 b);

int sign(float x);
