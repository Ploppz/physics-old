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


/* files */
#include "tmp.h"
#include "Renderer.h"
#include "shaders.h"
#include "glutils.h"
#include "Polygon.h"
#include "Geometry.h"
#include "BodySystem.h"
#include "World.h"
// Typewriter
#include "typewriter/FontTexture.h"
#include "typewriter/FontRenderer.h"

void error_callback(int error, const char* description);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

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

float zoom = 1;

bool left_down, right_down, up_down, down_down;

int main()
{
    // srand (time(NULL));
    // srand(10);
	GLFWwindow *window;
	GLFW_boilerPlate(&window, error_callback);
	glfwSetKeyCallback(window, key_callback);
    glfwSetScrollCallback(window, scroll_callback);

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
        p.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (300 + a), sin(-i*2.0f/numEdges * M_PI) * (300 + a)));
    }
	numEdges = 5;
	for (int i = 0; i < numEdges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * 100;
		q.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (60 + a), sin(-i*2.0f/numEdges * M_PI) * (60 + a)));
	} 
	for (int i = 0; i < numEdges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * 100;
		r.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * (a), sin(-i*2.0f/numEdges * M_PI) * (a)));
	} 

	World world;
    Body b1, b2, b3;

	b1 = world.bodies.addBody();
	b1.shape() = p;

    b2 = world.bodies.addBody();
    b2.shape() = q;
    b2.position_type() = ABSOLUTE;

    b3 = world.bodies.addBody();
    b3.shape() = r;
    b2.position_type() = ABSOLUTE;
	
	std::vector<float> polygonBuffer;
	b1.addToBuffer(polygonBuffer);
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
    renderer.setColor1(0.5f, 0.1f, 0.1f);
    renderer.setColor2(0.7f, 0.5f, 0.03f);

	int width, height;
	double ratio, timer;

	double mouseX, mouseY;
	
	while (!glfwWindowShouldClose(window))
	{
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

		// b1 approaches mouse pointer
		float xTarget = mouseX - width/2;
		float yTarget = - mouseY + height/2;
		// b1.velocity().x += (xTarget - b1.position().x) / 100;
		// b1.velocity().y += (yTarget - b1.position().y) / 100;
		// b1.velocity() *= 0.9;
        if (up_down) b3.position().y += 10;
        if (down_down) b3.position().y -= 10;
        if (right_down) b3.position().x += 10;
        if (left_down) b3.position().x -= 10;

		world.timestep(30);
// #if 0
        std::vector<Intersection> intersects = Polygon::overlaps(b2.shape(), b3.shape());
        if (intersects.size() > 0) {
            // renderer.setColor1(0.6f, 0.22f, 0.35f);
        } else {
            renderer.setColor1(0.5f, 0.1f, 0.1f);
        }
        for (Intersection i : intersects) {
            renderer.addDot(i.point);
        }
// #endif
        std::vector<Polygon> intersections = Polygon::intersection(b2.shape(), b3.shape());
        std::cout << "# " << intersections.size() << std::endl;
        for (auto it = intersections.begin(); it != intersections.end(); it ++) {
            // Body b = world.bodies.addBody();
            // b.shape() = *it;
            it->appendLinesToVector(renderer.lines_buffer);
        }
        // b1.shape().appendLinesToVector(renderer.lines_buffer);
        // b2.shape().appendLinesToVector(renderer.lines_buffer);
        // b3.shape().appendLinesToVector(renderer.lines_buffer);

        

		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		ratio = width / (float) height;
		timer = (float)glfwGetTime() * 2;

        renderer.render(width, height, zoom);


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

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    zoom *= (-yoffset / 10.f + 1);
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        switch(key)
        {
            case GLFW_KEY_UP:
            case GLFW_KEY_K:
                up_down = true;
                break;
            case GLFW_KEY_DOWN:
            case GLFW_KEY_J:
                down_down = true;
                break;
            case GLFW_KEY_LEFT:
            case GLFW_KEY_H:
                left_down = true;
                break;
            case GLFW_KEY_RIGHT:
            case GLFW_KEY_L:
                right_down = true;
                break;
        }
    }
    if (action == GLFW_RELEASE) {
        switch(key)
        {
            case GLFW_KEY_UP:
            case GLFW_KEY_K:
                up_down = false;
                break;
            case GLFW_KEY_DOWN:
            case GLFW_KEY_J:
                down_down = false;
                break;
            case GLFW_KEY_LEFT:
            case GLFW_KEY_H:
                left_down = false;
                break;
            case GLFW_KEY_RIGHT:
            case GLFW_KEY_L:
                right_down = false;
                break;
        }
    }
}
