#pragma once

#include <glm/glm.hpp>
#include <iostream>

#include "BodySystem.h"
#include "Polygon.h"

class Body;
struct EdgePoint;

std::ostream &operator << (std::ostream &lhs, const glm::vec2 &rhs);
std::ostream &operator << (std::ostream &lhs, const glm::vec3 &rhs);

float angle(glm::vec2 a, glm::vec2 b, glm::vec2 c);

float length_squared(glm::vec2 v);
bool zero_length(glm::vec2 v);
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

/* Distance from polygon */
float distance(glm::vec2 point, Polygon& p, int& out_closest_edge, float& out_closest_edge_alpha);
float distance(glm::vec2 point, Polygon& p);
bool distance_along_line(glm::vec2 point, glm::vec2 line_direction, Polygon& p, float& out_distance);
// Without polygon transformation
float distance_model(glm::vec2 point, Polygon& p, int& out_closest_edge, float& out_closest_edge_alpha);

float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end, float &out_alpha);
float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end);

float distance_line_segment(glm::vec2 p, glm::vec2 line_start, glm::vec2 line_end, float &out_alpha);

// Midpoint of two points
glm::vec2 middle(glm::vec2 a, glm::vec2 b);

/* Projection */
glm::vec2 project(glm::vec2 p, glm::vec2 line_a, glm::vec2 line_b);
float project(glm::vec2 point, glm::vec2 direction);


float signedArea(glm::vec2 a, glm::vec2 b, glm::vec2 c);

glm::vec3 barycentric(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 p);

// following 2: segments no lines..
bool intersect(glm::vec2 line1_a, glm::vec2 line1_b, glm::vec2 line2_a, glm::vec2 line2_b);
bool intersect(glm::vec2 line1_a, glm::vec2 line1_b, glm::vec2 line2_a, glm::vec2 line2_b, glm::vec2& point_of_intersection, float& alpha1, float& alpha2);
bool intersect_line_line_segment(glm::vec2 line_start, glm::vec2 line_direction, glm::vec2 segment_start, glm::vec2 segment_end,
                                float& out_line_alpha);
// return the x value of the intersection between a line and a constant y value
float intersect_horizontal(glm::vec2 line_start, glm::vec2 line_direction, float y_constant, float &alpha_out);
float intersect_vertical(glm::vec2 line_start, glm::vec2 line_direction, float x_constant, float &alpha_out);
// Note to self ^^^^ I do like the say we use line_start and line_direction, generalizing for lines & line segments

/* Intersect between segment and untransformed polygon */
// Note: only returns the first result found
bool intersect_segment_polygon_model(glm::vec2 line_start, glm::vec2 line_end, Polygon& p, EdgePoint& result);

bool inside(glm::vec2 point, Polygon& p);

/* Using fewer transformations forth and back, aims to provide numerical stability at very small differences/numbers */
bool inside_stable(glm::vec2 point, Polygon& p);

// Non-transformed Polygon
bool inside_model(glm::vec2 point, Polygon& p);


float cross(glm::vec2 a, glm::vec2 b);

int sign(float x);

// transformation which rotates to align new_x_axis with the x axis
glm::mat2 rotate_coor_system(glm::vec2 new_x_axis);

glm::vec2 transform(glm::vec2 position, glm::vec2 translation, float orientation);
glm::vec2 detransform(glm::vec2 position, glm::vec2 translation, float orientation);
