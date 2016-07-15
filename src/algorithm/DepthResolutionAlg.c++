#if 0
#include <glm/glm.hpp>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "DepthResolutionAlg.h"
#include "geometry/linestrip/LineStrip.h"
#include "../Body.h"
#include "../render/Renderer.h"

extern Renderer *g_renderer;
const double PI = 3.141592653589793;

void DepthResolutionAlg::init()
{
    iterations = 0;
}

void DepthResolutionAlg::iteration_start()
{
    biggest_depth = 0;
    ++ iterations;
}

/**** Thoughts
Maybe not a good idea, that the Polygon doesn't get updated position when
we change it via Body.
****/
void DepthResolutionAlg::treat(Body b1, Body b2, float delta_time)
{
    const bool VISUAL_DEBUG = true;
    std::cout << "** TREAT **" << std::endl << std::endl;
    /* just_plot_movement(b2, b1, 30, 30);  */
    std::vector<Intersection> intersections =
    Polygon::extract_intersections(   b1.shape(), b2.shape(), bool(b1.mode), bool(b2.mode)); 
    if (intersections.size() == 0) return;

    for (auto it = intersections.begin(); it != intersections.end(); it ++)
    {

#if 0 // auxilliary_lines not defined in this class..;
        if (VISUAL_DEBUG) { 
            glm::vec2 normal = it->find_normal_wrt(&b1.shape());
            LineStrip ls1 = it->cast_shadow_on(&b1.shape(), normal);
            LineStrip ls2 = it->cast_shadow_on(&b2.shape(), - normal);
            ls1.append_lines_to_vector(auxilliary_lines, 1, 1, 1);
            ls2.append_lines_to_vector(auxilliary_lines, 1, 1, 1); 
        }
#endif

        DepthContact contact = it->get_contact(&b1.shape(), &b2.shape());
        resolve(b1, b2, contact);

        if (VISUAL_DEBUG && contact.depth > 0)
            g_renderer->extra_line_buffer.add_vector(contact.subj_point.point_t(), contact.ref_point.point_t() - contact.subj_point.point_t());



        std::cout << "Contact depth: " << contact.depth << std::endl;
        if (contact.depth > biggest_depth)
            biggest_depth = contact.depth;
    }
    b1.update_polygon_state();
    b2.update_polygon_state();
}

/** Find out which body should be moved the most... then move them **/
void DepthResolutionAlg::resolve(Body b1, Body b2, DepthContact contact)
{
    
    // NOTE: normal should point out of reference, which is the edge owner
    /** let b1 be reference and b2 be subject **/
    if (&b1.shape() != contact.ref_point.parent)
        std::swap(b1, b2);
    
    /** **/
    // TODO: Should we also consider whethe the dot product is negative?
    float ref_part = fabs(glm::dot(b1.velocity(), contact.normal));
    float subj_part = fabs(glm::dot(b2.velocity(), contact.normal));
    float distribution;
    if (ref_part == 0 && subj_part == 0)
        distribution = 0.5f;
    else
        distribution = 2 * atan2(ref_part, subj_part) / PI; // fair number in [0, 1]

    /** MOVE **/
    b1.position() -= distribution * contact.depth * contact.normal * CORRECTION_STRENGTH;
    b2.position() += (1 - distribution) * contact.depth * contact.normal * CORRECTION_STRENGTH;
}

bool DepthResolutionAlg::done()
{
    return (biggest_depth < DEPTH_OK_THRESHOLD
            || iterations > MAX_ITERATIONS);
}
#endif
