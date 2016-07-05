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

#include <fenv.h>
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
#include "shaders.h"
#include "glutils.h"
#include "Body.h"
#include "BodySystem.h"
#include "World.h"
#include "Input.h"
#include "error_handling.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"
#include "render/Renderer.h"
#include "geometry/Intersection.h"
#include "geometry/Polygon.h"
#include "geometry/geometry.h"
#include "debug/StatisticsCollection.h"
#include "debug/debug.h"
// Typewriter
#include "typewriter/FontTexture.h"
#include "typewriter/FontRenderer.h"

void error_callback(int error, const char* description);

Polygon create_polygon(int num_edges, float inner_size, float outer_size);

void set_up_test1(World& world, Body& to_be_controlled);

template <typename T>
std::ostream& operator<< (std::ostream&, std::vector<T>);

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
/* For debugging */
FontRenderer *g_font_renderer;
Renderer *g_renderer;
StatisticsCollection *g_statistics;


int main()
{
    DebugBeginC(false);
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	// Input::Init();

    /* GLFW  & Input */
	GLFWwindow *window;
	GLFW_boilerPlate(&window, error_callback);
	glfwSetKeyCallback(window, Input::key_callback);
	glfwSetScrollCallback(window, Input::scroll_callback);
	glfwSetMouseButtonCallback(window, Input::mouse_button_callback);

    /* ImGui */
    ImGui_ImplGlfwGL3_Init(window, false);
    ImGuiIO& io = ImGui::GetIO();
    ImGuiStyle& style = ImGui::GetStyle();
	style.FramePadding.y = 1;
	style.ItemSpacing.y = 2;


	srand(1013);


	World world;

    Body big;
    set_up_test1(world, big);

    StatisticsCollection statistics;
    g_statistics = &statistics;

	Renderer renderer(world.bodies);
	g_renderer = &renderer;
	g_font_renderer = renderer.get_font_renderer();
    renderer.set_render_flag(POLYGON_SHOW_VELOCITY);
    // renderer.set_render_flag(POLYGON_SHOW_VERTEX_NUMBERS); 
	renderer.set_color_1(0.1f, 0.1f, 0);
	renderer.set_color_2(0.34f, 0.3f, 0.3f);

	int width, height;
	float zoom = 1;
	float center_x = 0, center_y = 0;
	double ratio, timer;

	double mouseX, mouseY;
	
	/** CONFIG **/
	const bool INTERACTIVE_FRAME = true;
	const int FRAME_DURATION_MS = 0;
	const float DELTA_TIME = 0.6f; // kinda milliseconds / 10..
	const float acceleration = 12; 
    glfwSwapInterval(1);
	/** **/

	int space_counter = 0;
	while (!glfwWindowShouldClose(window))
	{
        auto start_time = std::chrono::high_resolution_clock::now();
		/** Handle input **/
        // bounding_box.velocity() = glm::vec2(0);
        ImGui_ImplGlfwGL3_NewFrame();
        if (!io.WantCaptureMouse) {
            Input::UpdateMouse(window);
        }

		zoom *= (-Input::scroll / 10.f + 1);
		Input::scroll = 0;
		center_x -= Input::mouse_drag_x * zoom;
		center_y -= Input::mouse_drag_y * zoom;
		/** **/
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);
		glfwGetCursorPos(window, &mouseX, &mouseY);

		if (Input::keys[GLFW_KEY_N])	big.velocity().x -= acceleration/10;
		if (Input::keys[GLFW_KEY_E])	big.velocity().y -= acceleration/10;
		if (Input::keys[GLFW_KEY_I])	big.velocity().y += acceleration/10;
		if (Input::keys[GLFW_KEY_O])	big.velocity().x += acceleration/10;

		/* renderer.write_distances_to(big.shape(), bounding_box.shape()); */
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
		
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		ratio = width / (float) height;
		timer = (float)glfwGetTime() * 2;

		renderer.render(center_x, center_y, width, height, zoom);
        renderer.render(statistics);
        ImGui::Render();



		// Done drawing
		glfwSwapBuffers(window);
		glfwPollEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(FRAME_DURATION_MS));
        { /* Count FPS */
            auto delta_time = std::chrono::high_resolution_clock::now() - start_time;
            long long delta_time_micros = std::chrono::duration_cast<std::chrono::microseconds>(delta_time).count();
            statistics.add_value("FPS", 1000000.f / delta_time_micros);
        }
	}

    ImGui_ImplGlfwGL3_Shutdown();
	return 0;
}



void set_up_test1(World& world, Body& to_be_controlled)
{
	Body big, small, bounding_box;

	bounding_box = world.bodies.add_body();
	bounding_box.shape() = create_polygon(4, 1000, 0);
	bounding_box.position_type() = ABSOLUTE;
	bounding_box.rotation() = 0.03f;
	// bounding_box.rotation() = 0;
	// bounding_box.velocity() = glm::vec2(-4,  0);
	bounding_box.mode = POLYGON_OUTSIDE;

	big = world.bodies.add_body(bounding_box);
	big.shape() = create_polygon(10, 100, 200);
	big.position_type() = RELATIVE;
	big.velocity() = glm::vec2(55, 0.02f);
	big.rotation() = -0.05f;
	big.rotation() = 0;

#define SMALL_POLYGON_EXISTS true
	if (SMALL_POLYGON_EXISTS) {
		small = world.bodies.add_body(bounding_box);
		small.shape() = create_polygon(10, 50, 150);
		small.position_type() = RELATIVE;
		small.velocity() = glm::vec2(-10, -0.02f);
		small.rotation() = 0.05f;
		small.rotation() = 0; 
		// small.position().x = 500;
	}

    // Other polygons
    Body other = world.bodies.add_body(bounding_box);
    other.shape() = create_polygon(10, 50, 150);
    other.position_type() = RELATIVE;
    /* other.position().y = 500; */

    other = world.bodies.add_body(bounding_box);
    other.shape() = create_polygon(10, 50, 150);
    other.position_type() = RELATIVE;
    other.velocity() = glm::vec2(1, 5);
    /* other.position().y = - 500;  */

    to_be_controlled = big;
}



Polygon create_polygon(int num_edges, float inner_size, float outer_size)
{
    Polygon p;
    float a;
	for (int i = 0; i < num_edges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * outer_size;
		p.vertices.push_back(glm::vec2(cos(-i*2.0f/num_edges * M_PI) * (inner_size + a), sin(-i*2.0f/num_edges * M_PI) * (inner_size + a)));
	} 
	p.calculate_shape_dependent_variables();
    return p;
}




/*******************/
void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

template <typename T>
std::ostream& operator<< (std::ostream& out, std::vector<T> vec)
{
    out << "vector(";
    for (T e : vec) {
        out << e << ", ";
    }
    out << ")";
    return out;
}

