#pragma once

#include <glm/glm.hpp>
#include <iostream>

#include "Polygon.h"
#include "BodySystem.h"

using namespace glm;
std::ostream &operator << (std::ostream &lhs, vec2 &rhs);
std::ostream &operator << (std::ostream &lhs, vec3 &rhs);

float angle(vec2 a, vec2 b, vec2 c);

// Returns true if a points left of b
bool leftof(vec2 a, vec2 b);

float distance(vec2 a, vec2 b);
// Distance from point to line
float distance(vec2 p, vec2 line_a, vec2 line_b);
// Midpoint of two points
vec2 middle(vec2 a, vec2 b);

vec2 project(vec2 p, vec2 line_a, vec2 line_b);

float signedArea(vec2 a, vec2 b, vec2 c);

// First 3 points are for triangle.
vec3 barycentric(vec2 a, vec2 b, vec2 c, vec2 p);


bool intersect(vec2 line1_a, vec2 line1_b, vec2 line2_a, vec2 line2_b);

bool inside(vec2 point, Body& body);
