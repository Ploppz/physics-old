/*
 *	TODO
 *	migrate overlap code to somewhere sensible and make it work there...
 *
 *  For now, make a function to split a polygon up in convex polygons,
 *  then maybe use SAT or something else to find the Manifolds.
 *  	Disable axes that are there only because of the splits.
 * 	  Later, maybe optimize by finding the collision points directly knowing the colliding edges?
 *
 *  Collide with boundaries of window
 */

/* std */
#include <iostream>
#include <cassert>
#include <vector>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <limits.h>
#include <math.h>


/* GL */
#include <GL/glew.h>
#include <GLFW/glfw3.h>

/* glm */
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>


/* files */
#include "tmp.h"
#include "Renderer.h"
#include "shaders.h"
#include "glutils.h"
#include "BodySystem.h"
#include "World.h"
#include "Input.h"
#include "geometry/Intersection.h"
#include "geometry/Polygon.h"
#include "geometry/geometry.h"
// Typewriter
#include "typewriter/FontTexture.h"
#include "typewriter/FontRenderer.h"

void error_callback(int error, const char* description);


template <typename T>
void printVector(std::vector<T> v)
{
	for (auto it = v.begin(); it != v.end(); it ++)
	{
		std::cout << *it << std::endl;
	}
}


// Fullwindow rectangle
GLfloat rectangle[] = {
	// x, y, tex_x, tex_y
	-1.0f, -1.0f, 0.0f, 1.0f,
	1.0f, -1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f, 0.0f,

	-1.0f, -1.0f,	0.0f, 1.0f, 
	1.0f, 1.0f, 1.0f, 0.0f,
	-1.0f, 1.0f, 0.0f, 0.0f
};


// TODO: Problems with bufferRef garbage. Maybe this is the problem? Make default constructor etc
FontRenderer *fontRenderer;


int main()
{
    // srand (time(NULL));
    // Input::Init();

	GLFWwindow *window;
	GLFW_boilerPlate(&window, error_callback);
	glfwSetKeyCallback(window, Input::key_callback);
    glfwSetScrollCallback(window, Input::scroll_callback);
    glfwSetMouseButtonCallback(window, Input::mouse_button_callback);


	// 
    Polygon p, q, r;
	int numEdges = 4;
	float a;
    for (int i = 0; i < numEdges; i ++) {
        float size = 1000;
        p.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * size, sin(-i*2.0f/numEdges * M_PI) * size));
    }
	numEdges = 10;
	for (int i = 0; i < numEdges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * 100;
        float size = 200;
		q.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (size + a), sin(-i*2.0f/numEdges * M_PI) * (size + a)));
	} 
	for (int i = 0; i < numEdges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * 100;
        float size = 100;
		r.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (size + a), sin(-i*2.0f/numEdges * M_PI) * (size/3 + a)));
	} 
    p.calculate_shape_dependent_variables();
    q.calculate_shape_dependent_variables();
    r.calculate_shape_dependent_variables();

	World world;
    Body big, small, bounding_box;

    bounding_box = world.bodies.add_body();
    bounding_box.shape() = p;
    bounding_box.position_type() = ABSOLUTE;
    bounding_box.rotation() = 0.03f;
    // bounding_box.rotation() = 0;
    // bounding_box.velocity() = glm::vec2(-4,  0);
    bounding_box.mode = POLYGON_OUTSIDE;

    big = world.bodies.add_body(bounding_box);
    big.shape() = q;
    big.position_type() = RELATIVE;
    big.velocity() = glm::vec2(10, 0);
    big.rotation() = -0.05f;
    big.rotation() = 0;

#define SMALL_POLYGON_EXISTS false
    if (SMALL_POLYGON_EXISTS) {
        small = world.bodies.add_body(bounding_box);
        small.shape() = r;
        small.position_type() = RELATIVE;
        small.velocity() = glm::vec2(-10, 0);
        small.rotation() = 0.05f;
        small.rotation() = 0; 
        small.position().x = 500;
    }
	
    Renderer renderer(world.bodies);
    renderer.set_render_flag(POLYGON_SHOW_VELOCITY);
    renderer.set_render_flag(POLYGON_SHOW_VERTEX_NUMBERS);
    renderer.set_color_1(0.1f, 0.1f, 0);
    renderer.set_color_2(0, 0.5f, 0);
    world.bodies.renderer = &renderer;

	int width, height;
    float zoom = 1;
    float center_x = 0, center_y = 0;
	double ratio, timer;

	double mouseX, mouseY;
	
    /** CONFIG **/
    const bool INTERACTIVE_FRAME = true;
    const int FRAME_DURATION_MS = 20;
    const float DELTA_TIME = 0.6f; // kinda milliseconds / 10..
    const float speed = 12; // this is in fact acceleration
    /** **/

    int space_counter = 0;
	while (!glfwWindowShouldClose(window))
	{
        /** Handle input **/

        Input::UpdateMouse(window);

        zoom *= (-Input::scroll / 10.f + 1);
        Input::scroll = 0;
        center_x -= Input::mouse_drag_x * zoom;
        center_y -= Input::mouse_drag_y * zoom;
        /** **/
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);
		glfwGetCursorPos(window, &mouseX, &mouseY);

        Body& a = big;
        if (SMALL_POLYGON_EXISTS) {
            Body& b = small; 
            if (Input::keys[GLFW_KEY_UP]) b.position().y += speed;
            if (Input::keys[GLFW_KEY_DOWN]) b.position().y -= speed; 
            if (Input::keys[GLFW_KEY_RIGHT]) b.position().x += speed; 
            if (Input::keys[GLFW_KEY_LEFT]) b.position().x -= speed;  
        }
        //
        if (Input::keys[GLFW_KEY_N])    a.velocity().x -= speed/10;
        if (Input::keys[GLFW_KEY_E])    a.velocity().y -= speed/10;
        if (Input::keys[GLFW_KEY_I])    a.velocity().y += speed/10;
        if (Input::keys[GLFW_KEY_O])    a.velocity().x += speed/10;

        if (INTERACTIVE_FRAME) {
            if (Input::keys[GLFW_KEY_SPACE]) {
                if (space_counter == 0 || space_counter > 10) {
                    world.timestep(DELTA_TIME);
                }
                space_counter ++;
            } else {
                space_counter = 0;
            }
        }
        /* a.velocity() *= 0.95f; */

        if ( ! INTERACTIVE_FRAME ) {
            world.timestep(DELTA_TIME);
        }
        /* std::vector<Intersection> intersections = Polygon::extract_intersections(bounding_box.shape(), big.shape(), false, false);
        if (intersections.size() > 0)
            renderer.append_lines_to_vector(intersections[0]); */

        

		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		ratio = width / (float) height;
		timer = (float)glfwGetTime() * 2;

        renderer.render(center_x, center_y, width, height, zoom);



		// Done drawing
		glfwSwapBuffers(window);
		glfwPollEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(FRAME_DURATION_MS));
	}

	return 0;
}


void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

