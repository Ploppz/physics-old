#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <list>

#include "geometry/Polygon.h"
#include "glutils.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"

#include "algorithm/DepthResolutionAlg.h"
#include "algorithm/TimeResolutionAlg.h"

class Body;


enum PositionType { ABSOLUTE, RELATIVE };
enum PolygonMode { POLYGON_INSIDE = 0, POLYGON_OUTSIDE = 1 };


class BodySystem
{
	friend class Body;
 /** METHODS **/
 public:
	BodySystem();
	void timestep(float delta);


	Body add_body(); // Returns index of new body
	Body add_body(Body parent);
	int num_bodies() const { return count; }
	Body get_body(int index);
	Polygon& get_polygon(int body_index);


	// 
	glm::mat3 construct_matrix(int index);

	std::vector<Intersect> overlaps(Body a, Body b);

 private:
	/*** Physical and Geometric Treatment ***/
	void treat_body_tree(Body root, float delta_time);


 public:
	// Just for drawing...
	void add_vector(glm::vec2 point, glm::vec2 vec);
	// Tests
	void visualize_shadows(Body, Body, std::vector<Intersection>&);
	void just_plot_movement(Body b1, Body b2, float total_time, int samples);

 /** MEMBERS **/
 public:
	float simulation_speed = 1;

	std::vector<float> auxilliary_lines;
 private:
	bool alternator = false;
	TimeResolutionAlg resolution_alg;
	int count;
	int flags = 0;


	// Body data
	std::vector<float> mass;

	std::vector<PositionType> position_type;

	std::vector<glm::vec2> position;
	std::vector<glm::vec2> past_position;
	std::vector<glm::vec2> velocity;
	std::vector<glm::vec2> force;

	std::vector<float> orientation;
    std::vector<float> past_orientation;
	std::vector<float> rotation;
	std::vector<float> torque;

	// Temporary solution: only one polygon each Body
	std::vector<Polygon> shape;

	// Hierarchical tree of bodies
	std::vector<Body> parent;
	std::vector<std::vector<Body>> children;

	std::vector<Body> top_level_bodies;
};
