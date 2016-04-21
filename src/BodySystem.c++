/** Thoughts
update_body could be a member function of Body?
should BS work with Body or int?
**/
#include "Renderer.h"
#include <cmath>
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <glm/gtx/string_cast.hpp>
#include <algorithm>

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
        
        update_body(i, delta_time);
        
        /* shape[i].position = position[i];
        shape[i].orientation = orientation[i]; */
	}

    for (Body p : top_level_bodies)
    {
        treat_body_tree(p, delta_time);
    }
}
void BodySystem::update_body(int index, float delta_time)
{
    position[index] += delta_time * velocity[index] * simulation_speed;
    orientation[index] += delta_time * rotation[index] * simulation_speed; 

    /* Copy transformation information for the Polygon object */
    shape[index].position = position[index];
    shape[index].orientation = orientation[index];
}


/** Treating bodies: penetration and physical **/

void BodySystem::treat_body_tree(Body root, float delta_time)
{
    // Children can't go outside of parent polygon:
    root.mode = POLYGON_OUTSIDE;
    for (int i = 0; i < children[root.index].size(); i ++)
    {
        Body child = children[root.index][i];
        child.mode = POLYGON_INSIDE;

        /* Test against parent */
        treat(child, root, delta_time);

        /* Test against siblings */
        for (int j = i+1; j < children[root.index].size(); j ++)
        {
            Body second_child = children[root.index][j];
            second_child.mode = POLYGON_INSIDE;

            treat(child, second_child, delta_time);
        }
        /* Recurse the tree */
        treat_body_tree(child, delta_time);
    }
}

bool reacted = false;
void BodySystem::treat(Body b1, Body b2, float delta_time)
{
    const bool VISUAL_DEBUG = false;
    /* just_plot_movement(b2, b1, 30, 30);  */
    std::vector<Intersection> intersections = Polygon::extract_intersections(   b1.shape(), b2.shape(),
                                                                                bool(b1.mode), bool(b2.mode)); 
    if (intersections.size() == 0) {
        std::cout << "YES" << std::endl;
        return;
    } else {
        std::cout << "NO" << std::endl;
    }
    Contact contact = find_earliest_contact(b1, b2, intersections, delta_time);
    std::cout << "\trewind time: " << contact.rewind_time << std::endl;
    if (contact.rewind_time == 0) {
        std::cout << "Ignoring collision" << std::endl;
        return;
    }


    if (VISUAL_DEBUG)
    {
        glm::vec2 p = contact.ref_point.point_t();
        renderer->add_dot(p);
        renderer->add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
    }

    /* if (reacted) return; */
    reacted = true;
    /* simulation_speed = 0;
    return;  */
    resolve_penetration(b1, b2, contact);
    physical_reaction(b1, b2, contact);
    if (separating_at(b1, b2, contact)) {
    
    } else {
    
    }
}
void BodySystem::resolve_penetration(Body b1, Body b2, Contact c)
{
    update_body(b1.index, -c.rewind_time);
    update_body(b2.index, -c.rewind_time);
}
void BodySystem::physical_reaction(Body A, Body B, Contact c)
{
    /* Maybe inconsistent: but we try to keep A the subject here */
    if (c.ref_point.parent == &A.shape()) {
        std::swap(A, B);
    }

    const float e = 0.1f; // coefficient of restitution
    glm::vec2 r_ortho_A, r_ortho_B;
    glm::vec2 v_AB = velocity_of_point(A, c.subj_point, r_ortho_A) - velocity_of_point(B, c.ref_point, r_ortho_B);

    float impulse = - (1 + e) * glm::dot(v_AB, c.normal) /
                    (  (1.f/A.shape().mass + 1.f/B.shape().mass) * glm::dot(c.normal, c.normal)
                        + std::pow(glm::dot(r_ortho_A, c.normal), 2) / A.shape().moment_of_inertia
                        + std::pow(glm::dot(r_ortho_B, c.normal), 2) / B.shape().moment_of_inertia  );
    A.apply_impulse(impulse * c.normal, c.subj_point);
    /* std::cout << "A mass " << A.shape().mass << std::endl;
    std::cout << "A inv mass " << A.shape().mass << std::endl;
    std::cout << "Impulse:  " << impulse << std::endl; */
    B.apply_impulse( - impulse * c.normal, c.ref_point);
}
inline glm::vec2 BodySystem::velocity_of_point(Body b, EdgePoint p, glm::vec2 &out_r_ortho)
{
    out_r_ortho = p.point_t() - b.position();
    out_r_ortho = glm::vec2(- out_r_ortho.y, out_r_ortho.x);
    return b.velocity() + b.rotation() * out_r_ortho;
}


bool BodySystem::separating_at(Body A, Body B, Contact c)
{
    if (c.ref_point.parent == &A.shape()) {
        std::swap(A, B);
    }
    glm::vec2 trash;
    glm::vec2 v_AB = velocity_of_point(A, c.subj_point, trash) - velocity_of_point(B, c.ref_point, trash);
}



Contact BodySystem::find_earliest_contact(Body b1, Body b2, std::vector<Intersection>& intersections, float delta_time)
{
    Contact earliest_contact;
    earliest_contact.rewind_time = 0;

    if (intersections.size() == 0) {
        return earliest_contact;
    }

    for (auto it = intersections.begin(); it != intersections.end(); it ++)
    {
        Contact c = calculate_contact(b1, b2, *it, delta_time);
        if (c.rewind_time > earliest_contact.rewind_time)
            earliest_contact = c;
    }
    return earliest_contact;
}

Contact BodySystem::calculate_contact(Body b1, Body b2, Intersection& intersection, float delta_time)
{
    assert(renderer);
    renderer->append_lines_to_vector(intersection);
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
                Contact current_contact = rewind_out_of(*it, b2, b1, delta_time);
                if (current_contact.rewind_time > best_contact.rewind_time) {
                    best_contact = current_contact;
                }
            } else if (it->owner == &b2.shape()) {
                Contact current_contact = rewind_out_of(*it, b1, b2, delta_time);
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
inline Contact BodySystem::rewind_out_of(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
    const bool VISUAL_DEBUG = false;
    assert(vertex.intersect == false);
    const float error_threshold = 0.1f;
    float error;
    float time_offset = delta_time / 2.f;
    float time_step_size =  delta_time / 2.f;
    int closest_edge = -1;
    float closest_edge_alpha = 0;
    /* if (will_separate_in_future(vertex, reference, subject, 2 * delta_time)) {
        std::cout << "It was found they will separate in the future" << std::endl;
        Contact null_contact;
        null_contact.rewind_time = 0;
        return null_contact;
    } */

    
    if (VISUAL_DEBUG)
    {
        glm::vec2 point_now = relative_pos(vertex.point, reference, subject, 0);
        renderer->add_dot(point_now);
        glm::vec2 point_then = relative_pos(vertex.point, reference, subject, - delta_time);
        renderer->add_dot(point_then);
    }
    int count = 0;
    const int max_count = 100;
    while (true)
    {
        count ++;
        if (count > max_count) break;
        glm::vec2 point_at_time = relative_pos(vertex.point, reference, subject, - time_offset);
        if (VISUAL_DEBUG)
            renderer->add_dot(point_at_time);

        error = distance(point_at_time, reference.shape(), closest_edge, closest_edge_alpha);
        /** PERFORMANCE? get velocity, create better estimate of new time offset **/

        if (reference.mode == POLYGON_OUTSIDE)
            error = - error; // outside is inside

        if (abs(error) <= error_threshold && error > 0) // TODO the check for > 0 isn't optimal.
            break;

        time_step_size *= 0.5f;
        if (error > 0)  time_offset -= time_step_size;
        else            time_offset += time_step_size;
    }
    if (error < 0) {
        std::cout << "They are still a bit inside" << std::endl;
    }
    Contact result;
    result.rewind_time = time_offset;
    result.normal = Polygon::Edge(closest_edge, &reference.shape()).normal_tr();
    result.ref_point = EdgePoint(closest_edge, closest_edge_alpha, &reference.shape());
    result.subj_point = EdgePoint(vertex.vertex, 0, &subject.shape());
    return result;
}

inline bool BodySystem::will_separate_in_future(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
    /** This function only considers the next frame, but to avoid critical situations,
        it should minimize intersection depth mid-frame.  **/
    glm::vec2 next_point = relative_pos(vertex.point, reference, subject, delta_time);
    float error = distance(next_point, reference.shape());
    return (error < 0);
}
void BodySystem::just_plot_movement(Body reference, Body subject, float total_time, int samples)
{
    Polygon& p = subject.shape();
    for (int vertex = 0; vertex < p.vertices.size(); vertex ++ )
    {
        for (int i = 0; i < samples; i ++)
        {
            float time = total_time / samples * i;
            glm::vec2 rel_pos = relative_pos(p.transformed(vertex), reference, subject, -time);
            renderer->add_dot(rel_pos);
        }
    }
}

inline glm::vec2 center(Body body, float time_offset)
{
    return body.position() + time_offset * body.velocity();
}

/** Notes on the math **/
/*
_Relative position_ is not atm necessary: the thought was that we can keep the reference still, but
it won't have any positive effect unless we pre-transform the reference polygon to its current state.
.. Which we will do: we will benefit from making a transformed copy of each incident polygon (?).
*/

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
    glm::vec2 relative_pos = translated_pos.x * unit_vector(- reference.rotation() * time_offset)
        + translated_pos.y * unit_vector_wrt_y(- reference.rotation() * time_offset)
        // constant (translation affiliated with rotation about a point):
        + reference.position() - unit_vector(- reference.rotation() * time_offset) * reference.position().x
        - unit_vector_wrt_y(- reference.rotation() * time_offset) * reference.position().y;
    // PERFORMANCE notice the repeated use of the angle ref.rotation() * time... ^
    return relative_pos;
}

/** BodySystem: Body management **/
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
    Body new_body = Body(this, count - 1);
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

//////////////
/*** Body ***/
//////////////

Body::Body(BodySystem *system, int index)
	: mode(POLYGON_INSIDE), system(system), index(index)
{ }

void Body::apply_impulse(glm::vec2 impulse, EdgePoint point)
{
    glm::vec2 r_ortho = point.point_t() - position();
    r_ortho = glm::vec2( - r_ortho.y, r_ortho.x);
    velocity() += impulse / shape().mass;
    rotation() += glm::dot(r_ortho, impulse) / shape().moment_of_inertia; 
}
void Body::apply_impulse(glm::vec2 impulse, glm::vec2 point)
{
    assert(!"Needs implemenation.");
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
