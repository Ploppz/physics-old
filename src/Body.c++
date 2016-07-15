#include <glm/glm.hpp>
/* src */
#include "Body.h"
#include "BodySystem.h"
#include "debug/debug.h"

using namespace glm;

//////////////
/*** Body ***/
//////////////

Body::Body(int index, BodySystem *system)
	: mode(POLYGON_INSIDE), system(system), index(index)
{ }

void Body::apply_impulse(glm::vec2 impulse, EdgePoint point)
{
	glm::vec2 r_ortho = point.point_t() - position();
	r_ortho = glm::vec2( - r_ortho.y, r_ortho.x);
	velocity() += impulse / shape().get_mass();
	rotation() += glm::dot(r_ortho, impulse) / shape().get_moment_of_inertia(); 
}
void Body::apply_impulse(glm::vec2 impulse, glm::vec2 point)
{
	assert(!"Needs implemenation.");
}
void Body::update(float delta_time)
{
	position() += delta_time * velocity();
	orientation() += delta_time * rotation(); 

	update_polygon_state();
}
void Body::update_polygon_state()
{
	shape().position = position();
	shape().orientation = orientation();
}
void Body::save_placement()
{
    ensure_valid();
    past_position() = position();
    past_orientation() = orientation();
}

vec2 Body::real_position()
{
	if (position_type() == ABSOLUTE || parent().index == -1) {
		return position();
	} else { // FIXED (relative to parent)
		return position() + parent().real_position();
	}
}

Body& Body::operator++ ()
{
    ++ index;
	return *this;
}
bool Body::is_valid()
{
	return index >= 0 && index < system->body_count;
}
void Body::ensure_valid()
{
	if (! is_valid()) {
		std::cerr << "Body is not valid: index: " << index << std::endl;
		exit(1);
	}
}

std::ostream& operator<< (std::ostream& out, Body b)
{
	out << "velocity\t=\t" << b.velocity()
		<< "\nposition\t=\t" << b.position()
		<< "\norientation\t=\t" << b.orientation()
		<< "\nrotation\t=\t" << b.rotation()
		<< std::endl;
	return out;
}

