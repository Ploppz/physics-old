#include <glm/glm.hpp>
#include <cmath>
#include <cfloat>
#include <algorithm>

/* src */
#include "debug/debug.h"
#include "geometry.h"
#include "EdgePoint.h"

using namespace glm;

float angle(vec2 a, vec2 b, vec2 c)
{
    a -= b;
    c -= b;
    return acos(dot(a, c) / (length(a) * length(c)));
}
float length_squared(vec2 v)
{
    return v.x*v.x + v.y*v.y;
}
bool zero_length(vec2 v)
{
    return v.x == 0 && v.y == 0;
}
bool leftof(vec2 a, vec2 b)
{
    // Cross product
    return (a.x * b.y - a.y * b.x) < 0;
}

float distance(vec2 a, vec2 b)
{
	vec2 diff = a - b;
	return sqrt(diff.x * diff.x + diff.y * diff.y);
}
float distance_line(vec2 p, vec2 line_a, vec2 line_b)
{
	return distance(p - line_a, project(p, line_a, line_b));
}
float distance_line_segment(vec2 p, vec2 line_start, vec2 line_end)
{
    float length_sq = length_squared(line_start - line_end);
    if (length_sq == 0.f) return distance(p, line_start);
    // clamp t between 0 and 1 because it's a _line segment_
    const float t = std::max(0.f, std::min(1.f, dot(p - line_start, line_end - line_start) / length_sq));
    const vec2 projection = line_start + t * (line_end - line_start);
    return distance(p, projection);
}
float distance_line_segment(vec2 p, vec2 line_start, vec2 line_end, float &out_alpha)
{
    float length_sq = length_squared(line_start - line_end);
    if (length_sq == 0.f) return distance(p, line_start);
    // clamp t between 0 and 1 because it's a _line segment_
    out_alpha = std::max(0.f, std::min(1.f, dot(p - line_start, line_end - line_start) / length_sq));
    const vec2 projection = line_start + out_alpha * (line_end - line_start);
    return distance(p, projection);
}
vec2 middle(vec2 a, vec2 b)
{
	return a + ((b - a) / 2.0f);
}
// Returns (projected p) - line_a?
vec2 project(vec2 p, vec2 line_a, vec2 line_b)
{
	vec2 u = line_b - line_a;
	p = p - line_a;
	return( dot(p, u)/dot(u, u) * u);
}

// Get projected size ("shadow") of polygon along direction
float project(Polygon polygon, vec2 direction)
{
    direction = normalize(direction);
    float min_shadow = FLT_MAX;
    float max_shadow = FLT_MIN;
    for (auto it = polygon.vertices.begin(); it != polygon.vertices.end(); it ++) {
        float shadow = project(*it, direction);
        if (shadow < min_shadow)
            min_shadow = shadow;
        if (shadow > max_shadow)
            max_shadow = shadow;
    }
    return max_shadow - min_shadow;
}
float project(vec2 point, vec2 direction)
{
    direction = normalize(direction);
    return (dot(point, direction) / dot(direction, direction));
}

float signedArea(vec2 a, vec2 b, vec2 c)
{
	return 0.5f*(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
}
vec3 barycentric(vec2 a, vec2 b, vec2 c, vec2 p)
{
	float A = signedArea(a, b, c);
	return vec3(	signedArea(b, c, p)/A,
					signedArea(c, a, p)/A,
					signedArea(a, b, p)/A	);
}

bool intersect(vec2 line1_a, vec2 line1_b, vec2 line2_a, vec2 line2_b, vec2& point_of_intersection, float& alpha1, float& alpha2)
{
    vec2 a_diff = line2_a - line1_a;
    vec2 line1_vec = line1_b - line1_a;
    vec2 line2_vec = line2_b - line2_a;
    float cross_vec = cross(line1_vec, line2_vec);
    if (cross_vec == 0) return false; // This means colinear (if also cross(a_diff, line2_vec) == 0) OR parallel (else)
    alpha1 = cross(a_diff, line2_vec) / cross_vec; //TODO wrong order??
    alpha2 = cross(a_diff, line1_vec) / cross_vec;
    if (alpha1 > 0 && alpha2 > 0 && alpha1 < 1 && alpha2 < 1) {
        point_of_intersection = line1_a + (line1_vec * alpha1);
        return true;
    }
    return false;
}
bool intersect(vec2 line1_a, vec2 line1_b, vec2 line2_a, vec2 line2_b)
{
    // Lazy
    vec2 a;
    float b, c;
    return intersect(line1_a, line1_b, line2_a, line2_b, a, b, c);
}
float intersect_horizontal(vec2 line_start, vec2 line_direction, float y_constant, float &alpha_out)
{
    // Returns x value of intersection
    assert(line_direction.y != 0);
    alpha_out = (y_constant - line_start.y)/line_direction.y;
    return (line_start.x + alpha_out * line_direction.x);
}
float intersect_vertical(vec2 line_start, vec2 line_direction, float x_constant, float &alpha_out)
{
    // Returns x value of intersection
    alpha_out = (x_constant - line_start.x)/line_direction.x;
    return (line_start.y + alpha_out * line_direction.y);
}

bool intersect_segment_polygon_model(glm::vec2 line_start, glm::vec2 line_end, Polygon& p, EdgePoint &result)
{
    float alpha1, alpha2;
    vec2 point_of_intersection;
    // PERFORMANCE can loop with less transformation calculations..
    Polygon::Edge start(0, &p);
    Polygon::Edge it(0, &p);
    do {
        bool intersecting = intersect(line_start, line_end, it.start(), it.end(), // in
                point_of_intersection, alpha1, alpha2); // out
        if (intersecting) {
            result.parent = &p;
            result.index = it.get_index();
            result.alpha = alpha2;
            return true;
        }
        ++ it;
    } while (it != start);

    /* No intersection */
    return false;
}

std::ostream &operator << (std::ostream &lhs, const vec2 &rhs)
{
	lhs << "(" << rhs.x << ", " << rhs.y << ")";
	return lhs;
}
std::ostream &operator << (std::ostream &lhs, const vec3 &rhs)
{
	lhs << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	return lhs;
}

int sign(float x)
{
    return (x < 0) ? -1 : 1;
}

bool inside(vec2 point, Polygon& p)
{
    LineSegment edge;
    int counter = 0;
    glm::vec2 delta;
    float t;
    for (int i = 0; i < p.num_edges(); i ++)
    {
        edge = p.get_edge(i);
        //TODO did some changes here
        edge.first = p.transform(edge.first);
        edge.second = p.transform(edge.second);
        edge.first -= point;
        edge.second -= point;
        if (sign(edge.first.y) != sign(edge.second.y)) { // The edge crosses the x-axis
            // Find out where it intersects the x-axis
            delta = edge.second - edge.first;
            t = - edge.first.y / delta.y;
            if ((edge.first + delta * t).x > 0) {
                counter ++;
            }
        }
    }
    return counter % 2 == 1;
}
bool inside_model(vec2 point, Polygon& p)
{
    LineSegment edge;
    int counter = 0;
    glm::vec2 delta;
    float t;
    for (int i = 0; i < p.num_edges(); i ++)
    {
        edge = p.get_edge(i);
        edge.first -= point;
        edge.second -= point;
        if (sign(edge.first.y) != sign(edge.second.y)) { // The edge crosses the x-axis
            // Find out where it intersects the x-axis
            delta = edge.second - edge.first;
            t = - edge.first.y / delta.y;
            if ((edge.first + delta * t).x > 0) {
                counter ++;
            }
        }
    }
    return counter % 2 == 1;
}

/**
 * CCW: Positive on the outside, negative on the inside
**/
float distance(vec2 point, Polygon& p, int& out_closest_edge, float& out_closest_edge_alpha)
{
    assert( ! (isnan(point.x) || isnan(point.y)));
    // PERFORMANCE incorporate the inside 
    
    // Experimental way to iterate
    Polygon::Edge start(0, &p);
    Polygon::Edge it(0, &p);
    float min_distance = FLT_MAX;
    do
    {
        // PERFORMANCE We transform twice per iteration ..
        float alpha_on_edge;
        float distance_from_edge = distance_line_segment(point, it.start_tr(), it.end_tr(), alpha_on_edge);
        if (distance_from_edge < min_distance) {
            min_distance = distance_from_edge;
            out_closest_edge = it.get_index();
            out_closest_edge_alpha = alpha_on_edge;
        }
        ++ it;
    } while (it != start);
    assert (min_distance != FLT_MAX);
    bool is_inside = inside(point, p) ^ p.CCW;
    if (is_inside)
        min_distance = - min_distance;
    return min_distance;
}

float distance(vec2 point, Polygon& p)
{
    // PERFORMANCE: specialize function rather than reuse
    int out_closest_edge;
    float out_closest_edge_alpha;
    return distance(point, p, out_closest_edge, out_closest_edge_alpha);
}

float distance_model(glm::vec2 point, Polygon& p, int& out_closest_edge, float& out_closest_edge_alpha)
{
    DebugBegin();
    assert( ! (isnan(point.x) || isnan(point.y)));
    // PERFORMANCE incorporate the inside 
    
    // Experimental way to iterate ( PERFORMANCE )
    Polygon::Edge start(0, &p);
    Polygon::Edge it(0, &p);
    float min_distance = FLT_MAX;
    do
    {
        float alpha_on_edge;
        float distance_from_edge = distance_line_segment(point, it.start(), it.end(), alpha_on_edge);
        if (alpha_on_edge != 1) { // Disregard, and let the same happen for the next edge, with alpha=0
            if (distance_from_edge < min_distance) {
                min_distance = distance_from_edge;
                out_closest_edge = it.get_index();
                out_closest_edge_alpha = alpha_on_edge;
            }
        }
        ++ it;
    } while (it != start);
    assert (min_distance != FLT_MAX);

    bool is_inside = inside_model(point, p) ^ p.CCW;
    if (is_inside)
        min_distance = - min_distance;
    return min_distance;
}

float cross(glm::vec2 a, glm::vec2 b)
{
    return (a.x * b.y - a.y * b.x);
}
mat2 rotate_coor_system(vec2 new_x_axis)
{
    float align_transform_d[4] = {new_x_axis.x, -new_x_axis.y, new_x_axis.y, new_x_axis.x}; //note: column major
    mat2 align_transform = make_mat2x2(align_transform_d);
    return align_transform;
}
glm::vec2 transform(glm::vec2 position, glm::vec2 translation, float orientation)
{
    float c = cos(orientation);
    float s = sin(orientation);
    return glm::vec2(   c * position.x - s * position.y + translation.x,
                        s * position.x + c * position.y + translation.y);
}
glm::vec2 detransform(glm::vec2 position, glm::vec2 translation, float orientation)
{
    position -= translation;
    float c = cos(orientation);
    float s = sin(orientation);
    return glm::vec2(   + c * position.x + s * position.y,
                        - s * position.x + c * position.y);
}
