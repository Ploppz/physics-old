#include <glm/glm.hpp>
/* std */
#include <list>
/* src */
#include "ResolutionAlg.h"
#include "config.h"
#include "Body.h"
#include "geometry/geometry.h"
#include "render/Graphics.h"
#include "debug/debug.h"

#define VISUAL_DEPTH_RESOLUTION true
#define VISUAL_BODY_PROGRESSION false

struct SpeedBackup {
    SpeedBackup(Body& b1, Body& b2) : b1(b1), b2(b2) {}
    void backup() {
        b1.velocity() = pre_velocity_b1; b2.velocity() = pre_velocity_b2;
        b1.rotation() = pre_rotation_b1; b2.rotation() = pre_rotation_b2;
    }
    Body& b1;
    Body& b2;
    glm::vec2 pre_velocity_b1 = b1.velocity();
    glm::vec2 pre_velocity_b2 = b2.velocity();
    float pre_rotation_b1 = b1.rotation();
    float pre_rotation_b2 = b2.rotation();
};

extern Graphics *g_graphics;
extern StatisticsCollection *g_statistics;

void ResolutionAlg::init()
{
    iterations = 0;
}
void ResolutionAlg::iteration_start()
{
    ++ iterations;
}
bool ResolutionAlg::done()
{
    return iterations >= MAX_ITERATIONS;
}

int count = 0;
void ResolutionAlg::treat(Body b1, Body b2, float delta_time)
{
    count ++;
    DebugBeginC(false);
    // dout << "Bodies " << b1.get_index() << " vs " << b2.get_index() << newl;
    std::vector<Intersection> intersections;
    intersections = Polygon::extract_intersections(b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 

    if (intersections.size() == 0) {
        dout << gray << "[No intersections.]" << newl;
        return;
    }

    SpeedBackup speed_backup(b1, b2);

    int time_iterations = 0;
    while (true) {
        { /* count */ 
            ++ time_iterations;
            if (time_iterations > MAX_TIME_ITER) {
                speed_backup.backup();
                break;
            }
        }

        TimeContact contact = find_earliest_contact_by_rewinding(b1, b2, intersections, delta_time);

        if (contact.rewind_time <= MIN_REWIND_TIME) {
            speed_backup.backup();
            break;
        }
        glm::vec2 rel_vel = relative_velocity(b1, b2, contact);
        if (glm::length(rel_vel) < MIN_REL_VEL) {
            speed_backup.backup();
            break;
        }

        dout << "Rewinding with rewind time = " << contact.rewind_time << newl; 
        _rewind(b1, b2, contact.rewind_time);
        physical_reaction(b1, b2, contact);
        if (separating_at(b1, b2, contact)) {
            dout << green << " -> Now Separating" << nocolor << newl;
            g_statistics->count("Now separating");
            unwind(b1, b2, contact.rewind_time);
        } else {
            dout << Red << " -> Now Approaching ********************************" << nocolor << newl;
            g_statistics->count("Now approaching");
            g_statistics->add_value("REL VEL", glm::length(rel_vel));
            unwind(b1, b2, contact.rewind_time);
            speed_backup.backup();
            break;
        }
        intersections = Polygon::extract_intersections(b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 
        if (intersections.size() == 0) {
            g_statistics->add_value("S time iterations", time_iterations);
            return;
        }

    }
    g_statistics->add_value("F time iterations", time_iterations);

    /* Rewinding failed -> revert and try again with Depth Resolution */
    /* Assumption is that at the previous frame, the polygons were *not overlapping* */
    /* And of course since it uses rewinding to check that, velocities/rotations must not be
     * changed earlier in the function.*/
    dout << " -> Rewinding failed... trying Depth Resolution" << newl << newl;
    g_statistics->count("Rewinding failed");
    treat_by_depth(b1, b2, delta_time);
}

void ResolutionAlg::treat_by_depth(Body b1, Body b2, float delta_time)
{
    DebugBeginC(false);

    int i = 0;

    std::list<DepthContact> contacts_resolved;
    std::vector<Intersection> intersections;
    intersections = Polygon::extract_intersections(b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 
    while (intersections.size() > 0){
        dout << "Iteration. Intersections = " << intersections.size() << " between bodies "
            << b1.get_index() << ", " << b2.get_index() << newl;

        // Loop through the vertices - find deepest contact //
        DepthContact deepest_contact;
        deepest_contact.depth = 0;
        for (Intersection& i : intersections)
        {
            for (HybridVertex& vertex : i.vertices)
            {
                if (!vertex.intersect) {
                    DepthContact contact = linear_find_contact(b1, b2, vertex);
                    if (contact.depth > deepest_contact.depth) {
                        deepest_contact = contact;
                    }
                }
            }
        }
        if (VISUAL_BODY_PROGRESSION) {
            g_graphics->append_lines(b1.shape());
            g_graphics->append_lines(b2.shape());
        }
        if (fabs(deepest_contact.depth) > INSIGNIFICANT_DEPTH) {
            resolve_by_depth(b1, b2, deepest_contact);
            contacts_resolved.push_back(deepest_contact); 
            if (VISUAL_BODY_PROGRESSION) {
                g_graphics->append_lines(b1.shape());
                g_graphics->append_lines(b2.shape());
            }
        } else {
            if (deepest_contact.depth == 0) {
                dout << Red << "Depth: Null contact returned. (" << deepest_contact.depth << ")" << newl;
                g_statistics->count("null contact returned");
            } else {
                g_statistics->count("depth_ignored");
            }
            break;
        }
        intersections = Polygon::extract_intersections(b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 

        ++ i;
        if (i >= MAX_DEPTH_ITER) {   
            break;
            dout.fatal("Max iterations.");
        }
    }
    g_statistics->add_value("depth iterations", i);

    for (DepthContact& contact : contacts_resolved)
    {
        physical_reaction(b1, b2, contact);
    }
}
DepthContact ResolutionAlg::linear_find_contact(Body subject, Body reference, HybridVertex& vertex)
{
    DebugBeginC(false);
    /* Ensure that subject is the owner of the vertex */
    if (&subject.shape() != vertex.owner) {
        std::swap(subject, reference);
    }

    /* Current point,          transformed to the coordinate system of _reference_ */
    glm::vec2 point_detransformed = reference.shape().detransform(vertex.point);
    /* Find point in the past, transformed to the coordinate system of _reference_ in the past */
    glm::vec2 past_point_detransformed = transform(subject.shape().model_vertex(vertex.vertex),
            subject.past_position(), subject.past_orientation());
    past_point_detransformed = detransform(past_point_detransformed,
            reference.past_position(), reference.past_orientation());

    /* Find intersect between the 'sweep line' and the reference Body */
    // NOTE/TODO: actually only need edge number of intersection. should it only return an int?
    EdgePoint intersect;
    bool intersect_happened =
        intersect_segment_polygon_model(past_point_detransformed, point_detransformed, reference.shape(), intersect);

    const bool EXTRA_INFO = false;
    if (EXTRA_INFO) {
        float dist;
        int closest_edge;
        float closest_edge_alpha;
        dist = distance_model(point_detransformed, reference.shape(), closest_edge, closest_edge_alpha);
        if (reference.mode == POLYGON_OUTSIDE) dist = - dist;
        dout << "Distance of current point: " << dist << newl;
        dist = distance_model(past_point_detransformed, reference.shape(), closest_edge, closest_edge_alpha);
        if (reference.mode == POLYGON_OUTSIDE) dist = - dist;
        dout << "Distance of past point: " << dist << newl;
    }
    if (!intersect_happened) {
        dout << red << "Intersect didn't happen - trying approx." << newl;
        return linear_find_contact_approx(subject, reference, vertex);
    }

    Polygon::Edge edge_of_intersect(intersect.index, intersect.parent);

    /* Find normal and depth by looking only at what edge is involved */
    DepthContact result;

    result.subj_point = EdgePoint(vertex.vertex, 0, &subject.shape());
    result.ref_point = intersect;

    // old way to obtain depth:
    result.depth = distance_line_segment(vertex.point, edge_of_intersect.start_tr(), edge_of_intersect.end_tr());
    result.normal = edge_of_intersect.normal_tr();
    result.ensure_normal_toward_subject();
    // return result;

    float depth;

    // Cast ray from subject point onto reference
    if ( distance_along_line(vertex.point, result.normal, reference.shape(), depth) ) {
        result.depth = depth;
    } else {
        // TODO Try the other approach here, e.g. normal * depth = distance between points.
        dout << red << "ray doesn't intersect polyogn?" << newl;
        g_statistics->count("ray failed");
        result.depth = 0;
        return result;
    }


    if (EXTRA_INFO) {
        dout << "Normal: " << result.normal << "   .. the direction in which Reference will move." << newl;
        dout << "Distance from reference to subject point: " <<
            (result.subj_point.point_t() - result.ref_point.point_t()) << newl;
    }

    dout << "End ... result.depth = " << result.depth << newl;
    return result;
}

DepthContact ResolutionAlg::linear_find_contact_approx(Body subject, Body reference, HybridVertex& vertex)
{
    DebugBeginC(false);
    /* Since this is called mostly in marginal situations..: */
    /* Just find the edge that is the closest to the point in the past */

    /* Current point,          transformed to the coordinate system of _reference_ */
    glm::vec2 point_detransformed = reference.shape().detransform(vertex.point);
    /* Find point in the past, transformed to the coordinate system of _reference_ in the past */
    glm::vec2 past_point_detransformed =
        transform(subject.shape().model_vertex(vertex.vertex), subject.past_position(), subject.past_orientation());
    past_point_detransformed =
        detransform(past_point_detransformed, reference.past_position(), reference.past_orientation());

    int closest_edge;
    float closest_edge_alpha;
    float dist;
    
    const bool EXTRA_INFO = false;
    { // test
        dist = distance_model(point_detransformed, reference.shape(), closest_edge, closest_edge_alpha);
        if (reference.mode == POLYGON_OUTSIDE) dist = - dist;
        if (EXTRA_INFO) dout << "Distance of current point: " << dist << newl;
    }
    dist = distance_model(past_point_detransformed, reference.shape(), closest_edge, closest_edge_alpha);
    if (reference.mode == POLYGON_OUTSIDE) dist = - dist;
    if (EXTRA_INFO) dout << "Distance of past point: " << dist << newl;


    Polygon::Edge edge_of_intersect(closest_edge, &reference.shape());
    EdgePoint edgepoint_of_intersect(closest_edge, closest_edge_alpha, &reference.shape());

    /* Get depth depending on the _current point_ */
    DepthContact result;
    result.subj_point = EdgePoint(vertex.vertex, 0, &subject.shape());
    result.ref_point = EdgePoint(closest_edge, closest_edge_alpha, &reference.shape());

    { // DEBUG
        dout << "subj point: " << vertex.vertex << " of " << subject.get_index() << newl;
        dout << "ref point: " << closest_edge << " of " << reference.get_index() << newl;
    }

    // Development notes:
    /* In non-approx method, depth was calculated by intersecting the normal ray with the subject polygon
     * --- I think that here, it's sufficient to do it the old way... TODO read this */
    result.depth =
        distance_line_segment(vertex.point, edge_of_intersect.start_tr(), edge_of_intersect.end_tr());

    /** Normal **/
    if (closest_edge_alpha == 0) {
        dout << Magenta << "Vertex vs. Vertex" << newl;
        /* Vertex vs. Vertex ---> normal is in direction of distance between points */
        result.normal = subject.shape().transformed(vertex.vertex) - edgepoint_of_intersect.point_t();
    } else {
        /* Else ---> normal is just the normal of the edge */
        result.normal = edge_of_intersect.normal_tr();
    }
    if (zero_length(result.normal)) {
        result.depth = 0;
        return result;
    } else {
        result.normal = glm::normalize(result.normal);
    }


    result.ensure_normal_toward_subject();
    /* Attempt: if subj_point is actually on the outside, flip the normal. */
    dist = distance(result.subj_point.point_t(), reference.shape()); 
    if (reference.mode == POLYGON_OUTSIDE) dist = - dist;
    if (dist >= 0) {
        // STATISTICS
        dout << green << "Flipped normal." << newl;
        result.normal = - result.normal; 
        // result.depth = 0;
    } else {
        dout << green << "Normal is ok." << newl;
    }

    if (EXTRA_INFO) dout << "Depth: " << result.depth << newl;
    if (EXTRA_INFO) dout << "Resulting point on reference: " << result.ref_point.point_t() << newl;
    return result;
}


const double PI = 3.141592653589793;

void ResolutionAlg::resolve_by_depth(Body subject, Body reference, DepthContact contact)
{
    DebugBeginC(false);
    /* Assumption: contact: normal points from ref_point to subj_point */
    /* Ensure that subject and contact.subj_point refer to the same polygon */
    if (&subject.shape() != contact.subj_point.parent)
        std::swap(subject, reference);
    const float strength = 1;
    const float additional_depth = 0.1f;
    dout << "Depth: " << contact.depth << newl;
    dout << "Subj point: " << contact.subj_point << newl;
    dout << "Ref point: " << contact.ref_point << newl;
    dout << "Resoution direction: " << contact.normal * contact.depth << newl;
    dout << " -- Subject = " << subject.get_index() << newl;

    /* LOG
    velocity direction + mass & just mass worked well on intended area 
      BUT... it makes several collisions work bad - it won't let smaller polygons 'claim their space'
      --> Could be fixed with a better multiple-collisions alg.?
    */

    /* Distribution */
    float distribution;
    { /* By velocity direction */
        float subj_part = fabs(glm::dot(subject.velocity(), contact.normal))  / subject.shape().get_mass();
        float ref_part = fabs(glm::dot(reference.velocity(), contact.normal)) / reference.shape().get_mass();
        if (subj_part == 0 && ref_part == 0)
            distribution = 0.5f;
        else
            distribution = 2 * atan2(subj_part, ref_part) / PI; // fair number in [0, 1]
    }
    // .. overrided by:
    { /* By mass */
        distribution = reference.shape().get_mass() / (subject.shape().get_mass() + reference.shape().get_mass());
    }

    subject.position() -=   (distribution)      * contact.normal * (contact.depth + additional_depth);
    reference.position() += (1 - distribution)  * contact.normal * (contact.depth + additional_depth);
    subject.update_polygon_state();
    reference.update_polygon_state();

    if (VISUAL_DEPTH_RESOLUTION)
    {
        g_graphics->set_line_color(glm::vec3(0.8f, 0, 0));
        g_graphics->line_renderer.add_vector(contact.subj_point.point_t(), - distribution * contact.normal * contact.depth); 
        g_graphics->line_renderer.add_dot(contact.subj_point.point_t(), distribution * contact.depth / 15.f); 

        g_graphics->set_line_color(glm::vec3(0.8f, 0.8f, 0));
        g_graphics->line_renderer.add_vector(contact.ref_point.point_t(), (1 - distribution) * contact.normal * contact.depth); 
        g_graphics->line_renderer.add_dot(contact.ref_point.point_t(), (1 - distribution) * contact.depth / 15.f);
    }
}



TimeContact ResolutionAlg::find_earliest_contact_by_rewinding(Body b1, Body b2, std::vector<Intersection>& intersections, float delta_time)
{
    DebugBegin();

	/** First check whether any one vertex cannot be rewound **/
	for (auto it = intersections.begin(); it != intersections.end(); it ++)
	{
		if (!rewindable(b1, b2, *it, delta_time)) {
			TimeContact null_contact;
			null_contact.rewind_time = 0;
            dout << "Null contact returned";
			return null_contact;
		}
	}
    dout << "Rewindable.";

	/** **/
	TimeContact earliest_contact;
	earliest_contact.rewind_time = 0;

	if (intersections.size() == 0) {
		return earliest_contact;
	}

	for (auto it = intersections.begin(); it != intersections.end(); it ++)
	{
		TimeContact c = calculate_contact_by_rewinding(b1, b2, *it, delta_time);
		if (c.rewind_time > earliest_contact.rewind_time)
			earliest_contact = c;
	}
	return earliest_contact;
}
bool ResolutionAlg::rewindable(Body b1, Body b2, Intersection& intersection, float delta_time)
{
    DebugBegin();
	for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
	{
		if (!separate_last_frame(*it, b1, b2, delta_time)) {
            dout << "Separate last frame = false, for some vertex";
			return false;
		}
	}
	return true;
}

TimeContact ResolutionAlg::calculate_contact_by_rewinding(Body b1, Body b2, Intersection& intersection, float delta_time)
{
    DebugBegin();
	assert(g_graphics);
	/* g_graphics->extra_line_buffer.append_lines_to_vector(intersection); */
	/** Possibilities
	 Send 'intersection' to rewind_out_of to only consider those edges?
	**/

	/* Find the TimeContact with the highest rewind time */
	TimeContact best_contact;
	best_contact.rewind_time = 0;
	for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
	{
		if ( ! (it->intersect)) {
			if (it->owner == &b1.shape()) {
				TimeContact current_contact = rewind_out_of(*it, b2, b1, delta_time);
				if (current_contact.rewind_time > best_contact.rewind_time) {
					best_contact = current_contact;
				}
			} else if (it->owner == &b2.shape()) {
				TimeContact current_contact = rewind_out_of(*it, b1, b2, delta_time);
				if (current_contact.rewind_time > best_contact.rewind_time) {
					best_contact = current_contact;
				}
			}
		}
	}
    if (best_contact.rewind_time == 0) dout << "Rewind time still 0" << newl;
	return best_contact;
}
inline TimeContact ResolutionAlg::rewind_out_of(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
    DebugBegin();
	assert(vertex.intersect == false);

	const bool VISUAL_DEBUG = false;

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
	if ( ! separate_last_frame(vertex, reference, subject, delta_time)) {
		dout << "It was found that they were not separate last frame" << newl;
		TimeContact null_contact;
		null_contact.rewind_time = 0;
		return null_contact;  
	}

	
	if (VISUAL_DEBUG)
	{
		glm::vec2 point_now = relative_pos(vertex.point, reference, subject, 0);
		g_graphics->line_renderer.add_dot(point_now);
		glm::vec2 point_then = relative_pos(vertex.point, reference, subject, - delta_time);
		g_graphics->line_renderer.add_dot(point_then);
	}
	int count = 0;
	while (true)
	{
		count ++;
		if (count > MAX_REWIND_ITER) break;
		glm::vec2 point_at_time = relative_pos(vertex.point, reference, subject, - time_offset);
		if (VISUAL_DEBUG)
			g_graphics->line_renderer.add_dot(point_at_time);

		error = distance(point_at_time, reference.shape(), closest_edge, closest_edge_alpha);
		/** PERFORMANCE? get velocity, create better estimate of new time offset **/

		if (reference.mode == POLYGON_OUTSIDE)
			error = - error; // outside is inside

		if (abs(error) <= ERROR_THRESHOLD && error > 0) // TODO the check for > 0 isn't optimal.
			break;

		time_step_size *= 0.5f;
		if (error > 0)  time_offset -= time_step_size;
		else			time_offset += time_step_size;
	}
	if (error < 0) {
		dout << "They are still a bit inside" << newl;
	}
	TimeContact result;
	result.rewind_time = time_offset;
	result.normal = Polygon::Edge(closest_edge, &reference.shape()).normal_tr();
	result.ref_point = EdgePoint(closest_edge, closest_edge_alpha, &reference.shape());
	result.subj_point = EdgePoint(vertex.vertex, 0, &subject.shape());
	return result;
}





void ResolutionAlg::_rewind(Body b1, Body b2, float time)
{
	b1.update(-time);
	b2.update(-time);
}
void ResolutionAlg::unwind(Body b1, Body b2, float time)
{
	b1.update(time);
	b2.update(time);
}
void ResolutionAlg::physical_reaction(Body A, Body B, TimeContact c)
{
	/* Maybe inconsistent: but we try to keep A the subject here */
	if (c.ref_point.parent == &A.shape()) {
		std::swap(A, B);
	}

	glm::vec2 r_ortho_A, r_ortho_B;
	glm::vec2 v_AB = velocity_of_point(A, c.subj_point, r_ortho_A) - velocity_of_point(B, c.ref_point, r_ortho_B);

	float impulse = - (1 + RESTITUTION) * glm::dot(v_AB, c.normal) /
					(  (1.f/A.shape().get_mass() + 1.f/B.shape().get_mass()) * glm::dot(c.normal, c.normal)
						+ std::pow(glm::dot(r_ortho_A, c.normal), 2) / A.shape().get_moment_of_inertia()
						+ std::pow(glm::dot(r_ortho_B, c.normal), 2) / B.shape().get_moment_of_inertia()  );
	A.apply_impulse(impulse * c.normal, c.subj_point);
	B.apply_impulse( - impulse * c.normal, c.ref_point);
}

glm::vec2 ResolutionAlg::relative_velocity(Body A, Body B, TimeContact c)
{
	glm::vec2 trash;
	return velocity_of_point(A, c.subj_point, trash) - velocity_of_point(B, c.ref_point, trash);
}
bool ResolutionAlg::separating_at(Body A, Body B, TimeContact c)
{
	if (c.ref_point.parent == &A.shape()) {
		std::swap(A, B); /* A is subject */
	}
	glm::vec2 v_AB = relative_velocity(A, B, c);
	/* Note to self: contact.normal should always be _out of_ reference */
	//TODO I'm here and the normal calculation takes into account CCW, which has been flipped and
	// it's now wrong.. question: does c.normal really reflect the _reference_?
	if ( (glm::dot(v_AB, c.normal) < 0) ^ B.mode) {
		return false;
	}
	return true;
}

inline glm::vec2 ResolutionAlg::velocity_of_point(Body b, EdgePoint p, glm::vec2 &out_r_ortho)
{
	out_r_ortho = p.point_t() - b.position();
	out_r_ortho = glm::vec2(- out_r_ortho.y, out_r_ortho.x);
	return b.velocity() + b.rotation() * out_r_ortho;
}


inline bool ResolutionAlg::will_separate_in_future(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
	/** This function only considers the next frame, but to avoid critical situations,
		it should minimize intersection depth mid-frame.  **/
	glm::vec2 next_point = relative_pos(vertex.point, reference, subject, delta_time);
	float error = distance(next_point, reference.shape());
	return (error < 0) ^ reference.mode ^ subject.mode;
}

// TODO PROBLEM is that subject is the bounding box
// and, uh, we don't consider which Body `vertex.point` belongs to
inline bool ResolutionAlg::separate_last_frame(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
    if (vertex.intersect) return true; // don't care
    // vertex should belong to subject.
    if (vertex.owner == &reference.shape())
        std::swap(reference, subject);
    DebugBeginC(false);
	glm::vec2 next_point = relative_pos(vertex.point, reference, subject, -delta_time);
	float error = distance(next_point, reference.shape());
    dout << "Error = " << error << newl;
    dout << "Modes = " << reference.mode << ", " << subject.mode << newl;
	/* std::cout << "Separation .. error  : " << error << std::endl;
	std::cout << "Separation .. mode   : " << reference.mode << ", " << subject.mode << std::endl; */
	return (error >= 0) ^ reference.mode ^ subject.mode;
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

glm::vec2 ResolutionAlg::relative_pos(glm::vec2 point, Body reference, Body subject, float time_offset)
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

