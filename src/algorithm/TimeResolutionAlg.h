#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <list>

#include "ResolutionAlg.h"
#include "Contact.h"
#include "geometry/EdgePoint.h"
#include "geometry/Intersection.h"

class Body;
/** Call tree..:
treat_body_tree
    treat
        find_earliest_contact
            calculate_contact
                rewind_out_of
                    relative_pos
        resolve_penetration
            update_body
        physical_reaction
            
**/

/* only used in a private method... also look in todo.md */
struct Placement {
    Placement(Body b, float rewind_time);
    glm::vec2 position;
    float orientation;
};

class TimeResolutionAlg
{
 public:
    void init();
    void iteration_start();
    void treat(Body, Body, float delta_time);
    bool done();
 public: /*private..*/
    void treat_by_depth(Body, Body, float delta_time);
    void resolve_by_depth(Body, Body, DepthContact contact);
    /* Uses previous frame to find incident edge & vertex, then simply calculates normal & depth */
    DepthContact linear_find_contact(Body subject, Body reference, HybridVertex& vertex,
                            Placement subj_past_placement, Placement ref_past_placement);


	void _rewind(Body b1, Body b2, float time);
    void unwind(Body b1, Body b2, float time);
	void physical_reaction(Body b1, Body b2, TimeContact c);
	bool separating_at(Body b1, Body b2, TimeContact c);

	TimeContact find_earliest_contact_by_rewinding(Body b1, Body b2, std::vector<Intersection>&, float time_since_last_update);
	TimeContact calculate_contact_by_rewinding(Body b1, Body b2, Intersection& b, float time_since_last_update);

	bool rewindable(Body b1, Body b2, Intersection& intersection, float time_since_last_update);

	inline bool will_separate_in_future(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_to_next_update);
	inline bool separate_last_frame(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_to_next_update);
	inline TimeContact rewind_out_of(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_since_last_update);
	static glm::vec2 relative_pos(glm::vec2 point, Body reference, Body subject, float time_offset);
	inline glm::vec2 velocity_of_point(Body b, EdgePoint p, glm::vec2 &out_r_ortho);


    // Unfinished, basic fallback in case Time method fails.
	std::list<TimeContact> simple_move_out_of(Body b1, Body b2, std::vector<Intersection>&);

 /** MEMBERS **/
 private:
    float biggest_depth;
    int iterations;
 private:
    const float CORRECTION_STRENGTH = 0.7f;
    const float DEPTH_OK_THRESHOLD = 10;
    const int MAX_ITERATIONS = 2;
};
