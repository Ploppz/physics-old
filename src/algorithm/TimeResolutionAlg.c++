#include "TimeResolutionAlg.h"
#include "../render/Renderer.h"
#include "Body.h"
#include <debug.h>


extern Renderer *g_renderer;

void TimeResolutionAlg::init()
{
    iterations = 0;
}
void TimeResolutionAlg::iteration_start()
{
    ++ iterations;
}
bool TimeResolutionAlg::done()
{
    return iterations >= MAX_ITERATIONS;
}
void TimeResolutionAlg::treat(Body b1, Body b2, float delta_time)
{
    DebugBegin();
    std::vector<Intersection> intersections =
            Polygon::extract_intersections(   b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 
    if (intersections.size() == 0){
        return;
    }

	TimeContact contact = find_earliest_contact_by_rewinding(b1, b2, intersections, delta_time);

	/* if (contact.rewind_time != 0) last_contact = contact;  */
	dout << "rewind time: " << contact.rewind_time << newl; 

	if (contact.rewind_time == 0) {
		dout << "SIMPLE DISPLACEMENT, WHICH DOESN'T WORK "<< newl;
		simple_move_out_of(b1, b2, intersections); 
		// physical reaction
	} else {
		resolve_penetration(b1, b2, contact);
		physical_reaction(b1, b2, contact);
		if (separating_at(b1, b2, contact)) {
			dout << "now Separating" << newl;
			b1.update(contact.rewind_time);
			b2.update(contact.rewind_time); 
			dout << "still Separating: " << std::boolalpha << separating_at(b1, b2, contact) << newl;
		} else {
			//
			dout << "***** now Approaching ********" << newl;
		}
	}
}


/** Private **/


TimeContact TimeResolutionAlg::find_earliest_contact_by_rewinding(Body b1, Body b2, std::vector<Intersection>& intersections, float delta_time)
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
bool TimeResolutionAlg::rewindable(Body b1, Body b2, Intersection& intersection, float delta_time)
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

TimeContact TimeResolutionAlg::calculate_contact_by_rewinding(Body b1, Body b2, Intersection& intersection, float delta_time)
{
    DebugBegin();
	assert(g_renderer);
	/* g_renderer->extra_line_buffer.append_lines_to_vector(intersection); */
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
inline TimeContact TimeResolutionAlg::rewind_out_of(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
    DebugBegin();
	assert(vertex.intersect == false);

	const bool VISUAL_DEBUG = false;
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
	if ( ! separate_last_frame(vertex, reference, subject, delta_time)) {
		dout << "It was found that they were not separate last frame" << newl;
		TimeContact null_contact;
		null_contact.rewind_time = 0;
		return null_contact;  
	}

	
	if (VISUAL_DEBUG)
	{
		glm::vec2 point_now = relative_pos(vertex.point, reference, subject, 0);
		g_renderer->extra_line_buffer.add_dot(point_now);
		glm::vec2 point_then = relative_pos(vertex.point, reference, subject, - delta_time);
		g_renderer->extra_line_buffer.add_dot(point_then);
	}
	int count = 0;
	const int max_count = 100;
	while (true)
	{
		count ++;
		if (count > max_count) break;
		glm::vec2 point_at_time = relative_pos(vertex.point, reference, subject, - time_offset);
		if (VISUAL_DEBUG)
			g_renderer->extra_line_buffer.add_dot(point_at_time);

		error = distance(point_at_time, reference.shape(), closest_edge, closest_edge_alpha);
		/** PERFORMANCE? get velocity, create better estimate of new time offset **/

		if (reference.mode == POLYGON_OUTSIDE)
			error = - error; // outside is inside

		if (abs(error) <= error_threshold && error > 0) // TODO the check for > 0 isn't optimal.
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





void TimeResolutionAlg::resolve_penetration(Body b1, Body b2, TimeContact c)
{
	b1.update(-c.rewind_time);
	b2.update(-c.rewind_time);
}
void TimeResolutionAlg::physical_reaction(Body A, Body B, TimeContact c)
{
	/* Maybe inconsistent: but we try to keep A the subject here */
	if (c.ref_point.parent == &A.shape()) {
		std::swap(A, B);
	}

	const float restitution = 0.5f; // coefficient of restitution
	glm::vec2 r_ortho_A, r_ortho_B;
	glm::vec2 v_AB = velocity_of_point(A, c.subj_point, r_ortho_A) - velocity_of_point(B, c.ref_point, r_ortho_B);

	float impulse = - (1 + restitution) * glm::dot(v_AB, c.normal) /
					(  (1.f/A.shape().mass + 1.f/B.shape().mass) * glm::dot(c.normal, c.normal)
						+ std::pow(glm::dot(r_ortho_A, c.normal), 2) / A.shape().moment_of_inertia
						+ std::pow(glm::dot(r_ortho_B, c.normal), 2) / B.shape().moment_of_inertia  );
	A.apply_impulse(impulse * c.normal, c.subj_point);
	B.apply_impulse( - impulse * c.normal, c.ref_point);
}
bool TimeResolutionAlg::separating_at(Body A, Body B, TimeContact c)
{
	if (c.ref_point.parent == &A.shape()) {
		std::swap(A, B); /* A is subject */
	}
	glm::vec2 trash;
	glm::vec2 v_AB = velocity_of_point(A, c.subj_point, trash) - velocity_of_point(B, c.ref_point, trash);
	/* Note to self: contact.normal should always be _out of_ reference */
	//TODO I'm here and the normal calculation takes into account CCW, which has been flipped and
	// it's now wrong.. question: does c.normal really reflect the _reference_?
	if ( (glm::dot(v_AB, c.normal) < 0) ^ B.mode) {
		return false;
	}
	return true;
}

inline glm::vec2 TimeResolutionAlg::velocity_of_point(Body b, EdgePoint p, glm::vec2 &out_r_ortho)
{
	out_r_ortho = p.point_t() - b.position();
	out_r_ortho = glm::vec2(- out_r_ortho.y, out_r_ortho.x);
	return b.velocity() + b.rotation() * out_r_ortho;
}

std::list<TimeContact> TimeResolutionAlg::simple_move_out_of(Body b1, Body b2, std::vector<Intersection>& intersections)
{
    DebugBegin();
	std::list<TimeContact> results;
	for (auto it = intersections.begin(); it != intersections.end(); it ++)
	{
		results.push_back(it->get_contact(  &b1.shape(), b1.mode == POLYGON_OUTSIDE,
											&b2.shape(), b2.mode == POLYGON_OUTSIDE));
	}
	int count = 0;
	do {
		if (count >= 10) {
			dout << "Reached max. iterations" << newl;
			break;
		}
		count ++;
		dout << "SimpleMoveOutOf iteration" << newl;

		Intersection& i = intersections[0];
		DepthContact contact = i.get_contact(  &b1.shape(), b1.mode == POLYGON_OUTSIDE,
											&b2.shape(), b2.mode == POLYGON_OUTSIDE);

		dout << "\tResulting depth: " << contact.depth << newl;
		/* The bodies are displaced weighed by their mass. contact.normal is assumed to be*/
		b1.position() -= contact.normal * (b2.shape().mass * contact.depth) / (b1.shape().mass + b2.shape().mass);
		b2.position() += contact.normal * (b1.shape().mass * contact.depth) / (b1.shape().mass + b2.shape().mass);

		b1.update_polygon_state();
		b2.update_polygon_state();

		dout << "=== Start extraction ... ===" << newl;
		intersections = Polygon::extract_intersections(   b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 
		dout << "Stop extraction ... " << newl;
	} while (intersections.size() > 0);
	dout << "SimpleMoveOutOf iterations: " << count << newl;
	return results;
}







/**
Numerically appoximate
at what time `point`, belonging to `subject`, no longer is inside `reference`
The Intersection is used as a guide of which edges of the reference to pick.
 -> In the future, we may want to check all edges or find a smarter way to select edges.
**/

inline bool TimeResolutionAlg::will_separate_in_future(HybridVertex vertex, Body reference, Body subject, float delta_time)
{
	/** This function only considers the next frame, but to avoid critical situations,
		it should minimize intersection depth mid-frame.  **/
	glm::vec2 next_point = relative_pos(vertex.point, reference, subject, delta_time);
	float error = distance(next_point, reference.shape());
	return (error < 0) ^ reference.mode ^ subject.mode;
}

// TODO PROBLEM is that subject is the bounding box
// and, uh, we don't consider which Body `vertex.point` belongs to
inline bool TimeResolutionAlg::separate_last_frame(HybridVertex vertex, Body reference, Body subject, float delta_time)
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

glm::vec2 TimeResolutionAlg::relative_pos(glm::vec2 point, Body reference, Body subject, float time_offset)
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
