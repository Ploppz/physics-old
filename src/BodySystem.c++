#include "Renderer.h"
#include <cmath>
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <glm/gtx/matrix_transform_2d.hpp>

#include "BodySystem.h"
#include "glutils.h"
#include "typewriter/FontRenderer.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"

extern FontRenderer *fontRenderer;

using namespace glm;

BodySystem::BodySystem()
	: count {}, mass {}, position {}, velocity {}, force {}, orientation {}, rotation {}, torque {}, shape {}
{
}

void BodySystem::timestep(float delta_time)
{
	for (int i = 0; i < count; i ++)
	{
		/* position[i] += delta_time * velocity[i];
        orientation[i] += delta_time * rotation[i]; */

        shape[i].position = position[i];
        shape[i].orientation = orientation[i];
	}
    for (Body p : top_level_bodies)
    {
        for (int i = 0; i < children[p.index].size(); i ++)
        {
            /* Test against parent */
            /* Test against siblings */
            for (int j = i; j < children[p.index].size(); j ++)
            {
            }
        }
    }
}
/* glm::mat3 BodySystem::construct_matrix(int index)
{
// TODO rotate
    return glm::translate(glm::mat3 {}, position[index]);
} */

Body BodySystem::add_body()
{
    position_type.push_back(ABSOLUTE);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
	rotation.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon() );

    children.push_back(std::vector<Body> {});
    parent.push_back(Body(this, -1)); // Invalid body

	++ count;
    Body new_body = Body(this,count - 1);
    top_level_bodies.push_back(new_body);
	return new_body;
}
Body BodySystem::add_body(Body parent)
{

    position_type.push_back(RELATIVE);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
	rotation.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon () );

    children[parent.index].push_back(Body(this, count));
    children.push_back(std::vector<Body> {});
    this->parent.push_back(parent);

	++ count;
	return Body(this, count - 1);
}

Body BodySystem::get_body(int index)
{
	return Body(this, index);
}

// Collision

typedef std::pair<glm::vec2, glm::vec2> LineSegment;

bool inside(vec2 point, Polygon polygon);



// For now this is copy-pasta code from Polygon: we just have to add the Body's pos.


// class BODY

Body::Body(BodySystem *system, int index)
	: system(system), index(index)
{
	// std::cout << "Body .. " << index << std::endl;
}

vec2 Body::real_position()
{
    if (position_type() == ABSOLUTE || parent().index == -1) {
        return position();
    } else { // RELATIVE
        return position() + parent().real_position();
    }
}

Body& Body::operator++ ()
{
    index ++;
    return *this;
}
bool Body::is_valid()
{
    return index >= 0 && index < system->count;
}


Contact BodySystem::calculate_earliest_contact(Body b1, Body b2, float time_since_last_update)
{
    std::vector<Intersection> intersections = Polygon::extract_intersections(b1.shape(), b2.shape(), false, false);
    Contact earliest_contact;

    if (intersections.size() == 0) {
        earliest_contact.rewind_time = 0;
        return earliest_contact;
    }

    for (auto it = intersections.begin(); it != intersections.end(); it ++)
    {
        Contact c = calculate_contact(b1, b2, *it, time_since_last_update);
        if (c.rewind_time > earliest_contact.rewind_time)
            earliest_contact = c;
    }
    return earliest_contact;
}

Contact BodySystem::calculate_contact(Body b1, Body b2, Intersection& intersection, float time_since_last_update)
{
    /** Possibilities
     Send 'intersection' to rewind_out_of to only consider those edges?
    **/

    /* Find the Contact with the highest rewind time */
    Contact best_contact;
    best_contact.rewind_time = 0;
    for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
    {
        if ( ! (it->intersect)) {
            if (it->owner == &b1.shape()) {
                std::cout << "b1" << std::endl;
                Contact current_contact = rewind_out_of(*it, b2, b1, time_since_last_update);
                if (current_contact.rewind_time > best_contact.rewind_time) {
                    best_contact = current_contact;
                }
            } else if (it->owner == &b2.shape()) {
                std::cout << "b2" << std::endl;
                Contact current_contact = rewind_out_of(*it, b1, b2, time_since_last_update);
                if (current_contact.rewind_time > best_contact.rewind_time) {
                    best_contact = current_contact;
                }
            }
        }
    }
    return best_contact;
}


/**
Numerically appoximate
at what time `point`, belonging to `subject`, no longer is inside `reference`
The Intersection is used as a guide of which edges of the reference to pick.
 -> In the future, we may want to check all edges or find a smarter way to select edges.
**/
inline Contact BodySystem::rewind_out_of(HybridVertex vertex, Body reference, Body subject, float time_since_last_update)
{
    assert(vertex.intersect == false);
    const float error_threshold = 1;
    float error;
    float time_offset = time_since_last_update / 2.f;
    float time_step_size =  time_since_last_update / 2.f;
    int closest_edge = -1;
    float closest_edge_alpha = 0;

    int count = 0;
    while (true)
    {
        count ++;
        if (count > 20) break;
        glm::vec2 point_at_time = relative_pos(vertex.point, reference, subject, - time_offset);
        renderer->add_dot(point_at_time);

        int closest_edge;
        error = distance(point_at_time, reference.shape(), closest_edge, closest_edge_alpha);
        // TODO get velocity, create better estimate of new time offset
        
        time_step_size *= 0.5f;
        if (abs(error) <= error_threshold)
            break;
        /* std::cout << count << ": " << error << std::endl; */

        if (error > 0)  time_offset -= time_step_size;
        else            time_offset += time_step_size;
    }
    Contact result;
    result.rewind_time = time_offset;
    result.normal = Polygon::Edge(closest_edge, &reference.shape()).normal_tr();
    result.ref_point = EdgePoint(closest_edge, closest_edge_alpha, &reference.shape());
    result.subj_point = EdgePoint(vertex.vertex, 0, &subject.shape());
    return result;
}
void BodySystem::just_plot_movement(Body reference, Body subject, float total_time, int samples)
{
    Polygon& p = subject.shape();
    for (int vertex = 0; vertex < p.vertices.size(); vertex ++ )
    {
        for (int i = 0; i < samples; i ++)
        {
            float time = total_time / samples * i;
            glm::vec2 rel_pos = relative_pos(p.transformed(vertex), reference, subject, time);
            renderer->add_dot(rel_pos);
        }
    }
}

inline glm::vec2 center(Body body, float time_offset)
{
    return body.position() + time_offset * body.velocity();
}
inline glm::vec2 BodySystem::relative_pos(glm::vec2 point, Body reference, Body subject, float time_offset)
{
    // PERFORMANCE There is a way to greatly optimize this.. a lot only has to be calculated once
    float angle_wrt_subject = angle_of_vector(point - subject.position());
    float angle_wrt_reference = angle_of_vector(point - reference.position());
    float distance_from_subject = glm::length(point - subject.position());
    float distance_from_reference = glm::length(point - reference.position());

    // First, absolute position:
    glm::vec2 translated_pos = subject.position() + subject.velocity() * time_offset
        + unit_vector(angle_wrt_subject + time_offset * subject.rotation()) * distance_from_subject;
    // Add translational velocity of reference
    // Translate:
    translated_pos -= time_offset * reference.velocity();
    // Add rotational velocity of reference
    // Rotate:
    glm::vec2 relative_pos = translated_pos.x * unit_vector(reference.rotation() * time_offset)
        + translated_pos.y * unit_vector_wrt_y(reference.rotation() * time_offset)
        // constant (translation affiliated with rotation about a point):
        + reference.position() - unit_vector(reference.rotation() * time_offset) * reference.position().x
        - unit_vector_wrt_y(reference.rotation() * time_offset) * reference.position().y;
    // PERFORMANCE notice the repeated use of the angle ref.rotation() * time... ^
    return relative_pos;
}
