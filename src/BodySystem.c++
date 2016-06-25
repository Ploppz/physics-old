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

#include "debug.h"
#include "BodySystem.h"
#include "Body.h"
#include "render/Renderer.h"
#include "geometry/geometry.h"
#include "geometry/linestrip/LineStrip.h"
#include "geometry/linestrip/LineStripSeries.h"
#include "glutils.h"
#include "typewriter/FontRenderer.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"
#include "algorithm/DepthResolutionAlg.h"

extern FontRenderer *fontRenderer;
extern Renderer *g_renderer;

using namespace glm;

BodySystem::BodySystem()
	: count {}, mass {}, position {}, velocity {}, force {}, orientation {}, rotation {}, torque {}, shape {}
{
}

void BodySystem::timestep(float delta_time)
{
    DebugBegin();
    delta_time *= simulation_speed;

	alternator = ! alternator;
	if (alternator) {
		g_renderer->extra_line_buffer.clear_buffer();
		for (int i = 0; i < count; i ++)
		{
			Body(i, this).update(delta_time);
		}
	} else {
		resolution_alg.init();

		do {
        dout << "ITER" << newl;
        resolution_alg.iteration_start();
		for (Body p : top_level_bodies)
		{
			treat_body_tree(p, delta_time);
		}
		} while (!resolution_alg.done());
	}
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

void BodySystem::visualize_shadows(Body b1, Body b2, std::vector<Intersection>& intersections)
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
}

void BodySystem::just_plot_movement(Body reference, Body subject, float total_time, int samples)
{
    Polygon& p = subject.shape();
    for (int vertex = 0; vertex < p.vertices.size(); vertex ++ )
    {
        for (int i = 0; i < samples; i ++)
        {
            float time = total_time / samples * i;
            glm::vec2 rel_pos = TimeResolutionAlg::relative_pos(p.transformed(vertex), reference, subject, -time);
            g_renderer->extra_line_buffer.add_dot(rel_pos);
        }
    }
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
	parent.push_back(Body(-1, this)); // Invalid body

	++ count;
	Body new_body = Body(count - 1, this);
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

	children[parent.index].push_back(Body(count, this));
	children.push_back(std::vector<Body> {});
	this->parent.push_back(parent);

	++ count;
	return Body(count - 1, this);
}

Body BodySystem::get_body(int index)
{
	return Body(index, this);
}



void BodySystem::add_vector(glm::vec2 point, glm::vec2 vec)
{
	const int radius = 4;
	const float arrow_angle = 2.4f;

	float vec_angle = atan2(vec.y, vec.x);
	glm::vec2 a1 = glm::vec2(cos(vec_angle - arrow_angle) * radius, sin(vec_angle - arrow_angle) * radius);
	glm::vec2 a2 = glm::vec2(cos(vec_angle + arrow_angle) * radius, sin(vec_angle + arrow_angle) * radius);

	auxilliary_lines.push_back(point.x);
	auxilliary_lines.push_back(point.y);
	auxilliary_lines.push_back(point.x + vec.x);
	auxilliary_lines.push_back(point.y + vec.y);

	auxilliary_lines.push_back(point.x + vec.x);
	auxilliary_lines.push_back(point.y + vec.y);
	auxilliary_lines.push_back(point.x + vec.x + a1.x);
	auxilliary_lines.push_back(point.y + vec.y + a1.y);

	auxilliary_lines.push_back(point.x + vec.x);
	auxilliary_lines.push_back(point.y + vec.y);
	auxilliary_lines.push_back(point.x + vec.x + a2.x);
	auxilliary_lines.push_back(point.y + vec.y + a2.y);

}
