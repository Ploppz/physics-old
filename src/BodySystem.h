#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <list>

#include "glutils.h"
#include "geometry/Polygon.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"

#include "algorithm/ResolutionAlg.h"
#include "algorithm/SAP.h"
#include "algorithm/PairOrderer.h"

class Body;


enum PositionType { FIXED, ABSOLUTE };
enum PolygonMode { POLYGON_INSIDE = 0, POLYGON_OUTSIDE = 1 };


class BodySystem
{
	friend class Body;
 public: /** Interface **/
	BodySystem();
	void timestep(float delta);


	Body add_body(); // Returns index of new body
	Body add_body(Body parent);
	int num_bodies() const { return body_count; }
	Body get_body(int index);
	Polygon& get_polygon(int body_index);

 private: /** Algorithm **/
	void treat_body_tree(Body root, float delta_time);
    void update_broadphase_alg();


 public: /** Tests/helpers **/
	void visualize_shadows(Body, Body, std::vector<Intersection>&);
	void just_plot_movement(Body b1, Body b2, float total_time, int samples);

 public: /** MEMBERS **/
	float simulation_speed = 1;

 private:
	ResolutionAlg resolution_alg;
    SAP<int> broadphase_alg;
    PairOrderer pair_orderer;
	int body_count;
	int flags = 0;
    // For debug
	bool alternator = false;


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

    // For broadphase
    std::vector<int> broadphase_id;
	// Hierarchical tree of bodies
	std::vector<Body> top_level_bodies;
	std::vector<Body> parent;
	std::vector<std::vector<Body>> children;
};
