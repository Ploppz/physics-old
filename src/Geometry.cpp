#include <glm/glm.hpp>
#include <cmath>

#include "Geometry.h"

using namespace glm;

float angle(vec2 a, vec2 b, vec2 c)
{
    a -= b;
    c -= b;
    return acos(dot(a, c) / (length(a) * length(c)));
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
float distance(vec2 p, vec2 line_a, vec2 line_b)
{
	return distance(p - line_a, project(p, line_a, line_b));
}
vec2 middle(vec2 a, vec2 b)
{
	return a + ((b - a) / 2.0f);
}
vec2 project(vec2 p, vec2 line_a, vec2 line_b)
{
	vec2 u = line_b - line_a;
	p = p - line_a;
	return( dot(p, u)/dot(u, u) * u);
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

bool intersect(vec2 line1_a, vec2 line1_b, vec2 line2_a, vec2 line2_b)
{
	// barycentric coordinates of line1_b.
	vec3 bary = barycentric(line1_a, line2_a, line2_b, line1_b);
	return bary.s <= 0 && bary.t >= 0 && bary.p >= 0;
}

std::ostream &operator << (std::ostream &lhs, vec2 &rhs)
{
	lhs << "(" << rhs.x << ", " << rhs.y << ")";
	return lhs;
}
std::ostream &operator << (std::ostream &lhs, vec3 &rhs)
{
	lhs << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	return lhs;
}
