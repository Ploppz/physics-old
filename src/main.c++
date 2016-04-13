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

bool file_exists(const std::string &name)
{
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}   
}

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

	// Typewriter
	FontTexture ft{};
	std::string s;
	if (file_exists("font_atlas")) {
		ft.useExistingAtlas("font_atlas", "metadata_atlas");
	} else {
		ft.generateAtlas("fonts/peep-07x14.bdf");

		ft.writeAtlasToFile("font_atlas", "metadata_atlas");
	}
    fontRenderer = new FontRenderer(1, ft);
    fontRenderer->setup();

	// 
    Polygon p, q, r;
	int numEdges = 5;
	float a;
    for (int i = 0; i < numEdges; i ++) {
        a = rand() / static_cast<float>(INT_MAX) * 100;
        float size = 300;
        p.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (size + a), sin(-i*2.0f/numEdges * M_PI) * (size + a)));
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

	World world;
    Body b1, b2, b3;

	// b1 = world.bodies.addBody();
	// b1.shape() = p;

    b2 = world.bodies.add_body();
    b2.shape() = q;
    b2.position_type() = ABSOLUTE;

    b3 = world.bodies.add_body();
    b3.shape() = r;
    b2.position_type() = ABSOLUTE;
	
	std::vector<float> polygonBuffer;
	// b1.addToBuffer(polygonBuffer);
    // int b2_offset = polygonBuffer.size();
    // b2.addToBuffer(polygonBuffer);

	// for (int i = 0; i < edges.size(); i ++)
	// {
		// addToBuffer(edges[i], polygonBuffer);
	// }


	// // Upload rect
	// GLuint polygonBufferRef= uploadVertices(polygonBuffer.data(), polygonBuffer.size() * sizeof(float));
	// // Attribute pointers
	// GLuint VAO = createVertexArrayObject();
	// glBindVertexArray(VAO);
	// glBindBuffer(GL_ARRAY_BUFFER, polygonBufferRef);
	// setFormat("position 2f color 3f", shaderProgram);

    
    Renderer renderer(world.bodies);
    renderer.set_color_1(0.1f, 0.1f, 0);
    renderer.set_color_2(0, 0.5f, 0);

	int width, height;
    float zoom = 1;
    float center_x = 0, center_y = 0;
	double ratio, timer;

	double mouseX, mouseY;
	
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
		// Reupload vertices
		// glBindBuffer(GL_ARRAY_BUFFER, polygonBufferRef);
		// BufferWriter<float> buffer(polygonBuffer.size());
		// b1.addToBuffer(buffer);
        // b2.addToBuffer(buffer);
		// glUnmapBuffer(GL_ARRAY_BUFFER);
		//
		glfwGetCursorPos(window, &mouseX, &mouseY);

        const float speed = 8;
        Body& b = b3;
        Body& a = b2;
        if (Input::keys[GLFW_KEY_UP]) {
            b.position().y += speed;
            // a.position().y -= speed;
        }
        if (Input::keys[GLFW_KEY_DOWN]) {
            b.position().y -= speed; 
            // a.position().y += speed;
        }
        if (Input::keys[GLFW_KEY_RIGHT]) {
            b.position().x += speed; 
            // a.position().x -= speed;
        }
        if (Input::keys[GLFW_KEY_LEFT]) {
            b.position().x -= speed; 
            // a.position().x += speed;
        }

		world.timestep(30);
        /* std::vector<Intersect> intersects = Polygon::overlaps(b2.shape(), b3.shape());
        for (Intersect i : intersects) {
            renderer.addDot(i.point);
        }  */
        std::vector<Intersection> intersections = Polygon::extract_intersections(b3.shape(), b2.shape(), false);
        if (intersections.size() > 0) {
        }
        // for (auto it = intersections.begin(); it != intersections.end(); it ++) {
            // it->appendLinesToVector(renderer.lines_buffer);
            // renderer.addDot(it->centroid());
        // }
        // b1.shape().appendLinesToVector(renderer.lines_buffer);
        // b2.shape().appendLinesToVector(renderer.lines_buffer);
        // b3.shape().appendLinesToVector(renderer.lines_buffer);

        

		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		ratio = width / (float) height;
		timer = (float)glfwGetTime() * 2;

        renderer.render(center_x, center_y, width, height, zoom);


		fontRenderer->render(width, height);
        fontRenderer->clearBuffer();

		// Done drawing
		glfwSwapBuffers(window);
		glfwPollEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(30));
	}

	return 0;
}


void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

