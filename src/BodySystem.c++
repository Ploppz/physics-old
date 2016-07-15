/** Thoughts
update_body could be a member function of Body?
should BS work with Body or int?
**/
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
#include <list>

/* src */
#include "BodySystem.h"
#include "Body.h"
#include "render/Renderer.h"
#include "geometry/geometry.h"
#include "glutils.h"
#include "typewriter/FontRenderer.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"
#include "algorithm/SAP.h"
#include "algorithm/PairManager.h"
#include "debug/debug.h"

extern FontRenderer *fontRenderer;
extern Renderer *g_renderer;
extern StatisticsCollection *g_statistics;

using namespace glm;

BodySystem::BodySystem()
	: body_count {}, mass {}, position {}, velocity {}, force {}, orientation {}, rotation {}, torque {}, shape {}
{
}
void BodySystem::timestep(float delta_time)
{
    g_statistics->count("frames");
    DebugBegin();
    dout << Green << "FRAME" << newl;
    
    delta_time *= simulation_speed;

	alternator = ! alternator;
#define alternate 0
#if alternate
	if (alternator) {
#endif
        /* Step forward */
		g_renderer->extra_line_buffer.clear_buffer();
		for (int i = 0; i < body_count; i ++)
		{
            Body(i, this).save_placement();
			Body(i, this).update(delta_time);
		}
#if alternate
	} else {
#endif
        /* Detect, resolve, react */
		resolution_alg.init();

        int i = 0;
		do {
            i ++;
            dout << "Iteration" << newl;
            resolution_alg.iteration_start();
            update_broadphase_alg();

            int active_pairs = 0;
            int used_pairs = 0;
            for (Pair pair : broadphase_alg.pairs)
            {
                ++ active_pairs;
                // dout << "Active pair: " << pair.box1 << " , " << pair.box2 << newl;
                Body b1(broadphase_alg.get_box_user_data(pair.box1), this);
                Body b2(broadphase_alg.get_box_user_data(pair.box2), this);
                // Treat only if they are siblings or have a parent-child relationship
                if (b1.parent() == b2.parent() || b1.parent() == b2 || b2.parent() == b1)
                {
                    if (b1.parent() == b2)
                        b2.mode = POLYGON_OUTSIDE;
                    else if (b2.parent() == b1)
                        b1.mode = POLYGON_OUTSIDE;


                    ++ used_pairs;
                    // dout << "TREATING " << b1.get_index() << ", " << b2.get_index() << newl;
                    resolution_alg.treat(b1, b2, delta_time);
                }
            
            }
            dout << "Active/used pairs: " << active_pairs << "/" << used_pairs << newl;
		} while (!resolution_alg.done());
        g_statistics->add_value("global iterations", i);
#if alternate
	}
#endif
    /* just_plot_movement(Body(0, this), Body(1, this), delta_time * 15, 15); */

}


/** Treating bodies: penetration and physical **/

void BodySystem::treat_body_tree(Body root, float delta_time)
{
    DebugBegin();
	// Children can't go outside of parent polygon:
	root.mode = POLYGON_OUTSIDE;
	for (int i = 0; i < children[root.index].size(); i ++)
	{
		Body child = children[root.index][i];
		child.mode = POLYGON_INSIDE;

		/* Test against parent */
		resolution_alg.treat(child, root, delta_time);

		/* Test against siblings */
		for (int j = i+1; j < children[root.index].size(); j ++)
		{
			Body second_child = children[root.index][j];
			second_child.mode = POLYGON_INSIDE;

			resolution_alg.treat(child, second_child, delta_time);
		}
		/* Recurse the tree */
		treat_body_tree(child, delta_time);
	}
}
void BodySystem::update_broadphase_alg()
{
    for (int i = 0; i < body_count; i ++)
    {
        AABB<2> aabb = shape[i].calc_bounding_box();
        g_renderer->extra_line_buffer.set_color(0.8f, 0.8f, 0.8f);
        g_renderer->extra_line_buffer.add_aabb(aabb.min[0], aabb.max[0], aabb.min[1], aabb.max[1]);
        broadphase_alg.update_box(broadphase_id[i], aabb);
    }
}

/* void BodySystem::visualize_shadows(Body b1, Body b2, std::vector<Intersection>& intersections)
{
	for (auto it = intersections.begin(); it != intersections.end(); it ++)
	{
		std::cout << *it << std::endl;
		// it->append_lines_to_vector(auxilliary_lines);
		glm::vec2 normal = it->find_normal_wrt(&b1.shape());
		add_vector(glm::vec2(0), normal * 15.f);
		std::cout << "*** Cast shadow on b1 ***" << std::endl;
		LineStrip ls1 = it->cast_shadow_on(&b1.shape(), normal);
		std::cout << "*** Cast shadow on b2 ***" << std::endl;
		LineStrip ls2 = it->cast_shadow_on(&b2.shape(), - normal);
		ls1.append_lines_to_vector(auxilliary_lines, 1, 1, 1);
		ls2.append_lines_to_vector(auxilliary_lines, 1, 1, 1); 
	}
} */

void BodySystem::just_plot_movement(Body reference, Body subject, float total_time, int samples)
{
    Polygon& p = subject.shape();
    for (Vertex v : p.vertices())
    {
        for (int i = 0; i < samples; i ++)
        {
            float time = total_time / samples * i;
            glm::vec2 rel_pos = ResolutionAlg::relative_pos(v.point, reference, subject, -time);
            g_renderer->extra_line_buffer.add_dot(rel_pos);
        }
    }
} 

/** BodySystem: Body management **/
Body BodySystem::add_body()
{
    DebugBegin();
	++ body_count;

	position_type.push_back(ABSOLUTE);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
    past_position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
    past_orientation.push_back( 0 );
	rotation.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon() );


	children.push_back(std::vector<Body> {});
	parent.push_back(Body(-1, this)); // Invalid body
	Body new_body = Body(body_count - 1, this);
	top_level_bodies.push_back(new_body);

    broadphase_id.push_back(broadphase_alg.add_box(AABB<2>(0, 0, 0, 0), new_body.get_index()));
    dout << "ADD BODY " << new_body.get_index() << newl;

	return new_body;
}
Body BodySystem::add_body(Body parent)
{
    DebugBegin();
	++ body_count;

	position_type.push_back(FIXED);
	mass.push_back( 0 );

	position.push_back( vec2(0, 0) );
    past_position.push_back( vec2(0, 0) );
	velocity.push_back( vec2(0, 0) );
	force.push_back( vec2(0, 0) );

	orientation.push_back( 0 );
    past_orientation.push_back( 0 );
	rotation.push_back( 0 );
	torque.push_back( 0 );
	shape.push_back( Polygon () );

	children[parent.index].push_back(Body(body_count, this));
	children.push_back(std::vector<Body> {});
	this->parent.push_back(parent);

    Body new_body = Body(body_count - 1, this);
    broadphase_id.push_back(broadphase_alg.add_box(AABB<2>(0, 0, 0, 0), new_body.get_index()));

    dout << "ADD CHILD BODY " << new_body.get_index() << newl;
	return new_body;
}

Body BodySystem::get_body(int index)
{
	return Body(index, this);
}
