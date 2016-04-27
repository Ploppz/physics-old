#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include <list>

#include "geometry/Polygon.h"
#include "glutils.h"
#include "geometry/geometry.h"
#include "geometry/Intersection.h"

/** Contact
rewind_time is positive for 'back in time'
**/

struct Contact {
    float rewind_time;
    // Collision normal in world coordinates
    glm::vec2 normal;

    // Collision incident points along the polygons
    EdgePoint ref_point;
    EdgePoint subj_point;
    Contact() {};
    Contact(IntersectionContact other) : rewind_time(0), normal(other.normal), ref_point(other.ref_point), subj_point(other.subj_point) {};
};

enum PositionType { ABSOLUTE, RELATIVE };

class Body;

class BodySystem
{
 friend class Body;
 public:
	BodySystem();
	void timestep(float delta);


    // fun
    void just_plot_movement(Body b1, Body b2, float total_time, int samples);

	Body add_body(); // Returns index of new body
    Body add_body(Body parent);
    int num_bodies() const { return count; }
	Body get_body(int index);
	Polygon& get_polygon(int body_index);


    // 
    glm::mat3 construct_matrix(int index);

    std::vector<Intersect> overlaps(Body a, Body b);
    // doesn't need..:
    Renderer* renderer;
    float simulation_speed = 1;
    Contact last_contact; // for visualization

    std::vector<float> auxilliary_lines;
 private:
    //
	int count;

    // Body data
	std::vector<float> mass;

    std::vector<PositionType> position_type;

	std::vector<glm::vec2> position;
	std::vector<glm::vec2> velocity;
	std::vector<glm::vec2> force;

	std::vector<float> orientation;
	std::vector<float> rotation;
	std::vector<float> torque;

	// Temporary solution: only one polygon each Body
	std::vector<Polygon> shape;

    // Hierarchical tree of bodies
    std::vector<Body> parent;
    std::vector<std::vector<Body>> children;

    std::vector<Body> top_level_bodies;


    void update_body(int body_index, float delta_time);

    /*** Physical and Geometric Treatment ***/
    void treat_body_tree(Body root, float delta_time);
    void treat(Body b1, Body b2, float delta_time);
    /* Returns time since resolved collision. */
    float resolve(Body b1, Body b2, std::vector<Intersection>& intersections, float delta_time);

    std::list<Contact> simple_move_out_of(Body b1, Body b2, std::vector<Intersection>&);
    Contact simple_make_contact(Body b1, Body b2, Intersection& i);


    void resolve_penetration(Body b1, Body b2, Contact c);
    void physical_reaction(Body b1, Body b2, Contact c);
    bool separating_at(Body b1, Body b2, Contact c);

    Contact find_earliest_contact_by_rewinding(Body b1, Body b2, std::vector<Intersection>&, float time_since_last_update);
    Contact calculate_contact_by_rewinding(Body b1, Body b2, Intersection& b, float time_since_last_update);

    bool not_rewindable(Body b1, Body b2, Intersection& intersection, float time_since_last_update);

    inline bool will_separate_in_future(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_to_next_update);
    inline bool separate_last_frame(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_to_next_update);
    inline Contact rewind_out_of(HybridVertex non_intersection_vertex, Body reference, Body subject, float time_since_last_update);
    inline glm::vec2 relative_pos(glm::vec2 point, Body reference, Body subject, float time_offset);
    inline glm::vec2 velocity_of_point(Body b, EdgePoint p, glm::vec2 &out_r_ortho);

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
};


enum PolygonMode {
    POLYGON_INSIDE = 0, POLYGON_OUTSIDE = 1
};

/** Body -- Smartpointer **/

class Body
{
	friend class BodySystem;
 public:
        PolygonMode mode;

	// Should invalidate it.
        Body() {mode = POLYGON_INSIDE;}
    // General
        void apply_impulse(glm::vec2 impulse, EdgePoint point);
        void apply_impulse(glm::vec2 impulse, glm::vec2 point); //TODO implement? don't know if I need it

        void update_polygon_state();

	// Access to members
        glm::vec2 real_position(); 
		/* inline float& mass() { return system->mass[index]; } */
        inline PositionType& position_type() { return system->position_type[index]; }
		inline glm::vec2& position() { ensure_valid(); return system->position[index]; }
		inline glm::vec2& velocity() { ensure_valid(); return system->velocity[index]; }
		inline glm::vec2& force() { ensure_valid(); return system->force[index]; }

		inline float& orientation() { ensure_valid(); return system->orientation[index]; }
		inline float& rotation() { ensure_valid(); return system->rotation[index]; }
		inline float& torque() { ensure_valid(); return system->torque[index]; }

		inline Polygon& shape() { ensure_valid(); return system->shape[index]; }

        inline Body& parent() { ensure_valid(); return system->parent[index]; }

	// Iterator?
        Body& operator++ ();
        bool is_valid();
        int get_index() {return index; }


 private:
	Body(BodySystem *system, int index);
    void ensure_valid();

	BodySystem *system;
	int index;
};
std::ostream& operator<< (std::ostream& out, Body b);
