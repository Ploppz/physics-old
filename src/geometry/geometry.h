#pragma once

#include <glm/glm.hpp>
#include <iostream>

#include "../BodySystem.h"
#include "Polygon.h"

class Body;

std::ostream &operator << (std::ostream &lhs, glm::vec2 &rhs);
std::ostream &operator << (std::ostream &lhs, glm::vec3 &rhs);

float angle(glm::vec2 a, glm::vec2 b, glm::vec2 c);
// Returns true if a points left of b
bool leftof(glm::vec2 a, glm::vec2 b);

float distance(glm::vec2 a, glm::vec2 b);
// Distance from point to line
float distance(glm::vec2 p, glm::vec2 line_a, glm::vec2 line_b);

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
