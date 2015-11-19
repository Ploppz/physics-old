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

Polygon p, q;
Body b1, b2;

// TODO: Problems with bufferRef garbage. Maybe this is the problem? Make default constructor etc
FontRenderer *fontRenderer;

float scale = 1;
float limit = 0;
int main()
{
    // srand (time(NULL));
    // srand(10);
	GLFWwindow *window;
	GLFW_boilerPlate(&window, error_callback, key_callback);

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
    fontRenderer->addText("The fox..lalalalalaLALA", 0, 0, false);
    fontRenderer->setup();

	// 
	int numEdges = 14;
	float a;
    for (int i = 0; i < numEdges; i ++) {
        a = rand() / static_cast<float>(INT_MAX) * 400;
        p.vertices.push_back(glm::vec2(cos(i*2.0f/numEdges * M_PI) * a, sin(i*2.0f/numEdges * M_PI) * a));
    }
    

    


    // p.vertices.push_back(glm::vec2(-300, 100));
    // p.vertices.push_back(glm::vec2(200, 100));
    // p.vertices.push_back(glm::vec2(50, 0));
    // p.vertices.push_back(glm::vec2(300, -100));
    // p.vertices.push_back(glm::vec2(-200, -100));
    // p.vertices.push_back(glm::vec2(-50, 0));


	numEdges = 10;
	for (int i = 0; i < numEdges; i ++) {
		a = rand() / static_cast<float>(INT_MAX) * 100;
		q.vertices.push_back(glm::vec2(cos(-i*2.0f/numEdges * M_PI) * a, sin(-i*2.0f/numEdges * M_PI) * a));
	} 



    std::vector<Triangle> triangles;

	std::vector<LineSegment> edges = p.decompose(triangles);
	
	BodySystem bodies {};
	World world;

	b1 = world.bodies.addBody();
	b1.shape() = p;
	b2 = world.bodies.addBody();
	b2.shape() = q;
	
	std::vector<float> polygonBuffer;
	b1.addToBuffer(polygonBuffer);
	// b2.addToBuffer(polygonBuffer);

	for (int i = 0; i < edges.size(); i ++)
	{
		addToBuffer(edges[i], b1, polygonBuffer);
	}

	// Compile shader
	GLuint vertexShader, fragmentShader; // unused
	GLuint shaderProgram  = createShaderProgram(shaders::shaders_vert, shaders::shaders_frag, vertexShader, fragmentShader);

	GLuint uni_proj = glGetUniformLocation(shaderProgram, "proj");
	GLuint uni_view = glGetUniformLocation(shaderProgram, "view");
	GLuint uni_model = glGetUniformLocation(shaderProgram, "model");
	

	// Upload rect
	GLuint polygonBufferRef= uploadVertices(polygonBuffer.data(), polygonBuffer.size() * sizeof(float));
	// Attribute pointers
	GLuint VAO = createVertexArrayObject();
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, polygonBufferRef);
	setFormat("position 2f color 3f", shaderProgram);

    // Triangle buffer
    
    std::vector<float> triangle_buffer;
    for (Triangle t : triangles) {
        triangle_buffer.push_back(t.a.x);
        triangle_buffer.push_back(t.a.y);
        triangle_buffer.push_back(t.color.r);
        triangle_buffer.push_back(t.color.g);
        triangle_buffer.push_back(t.color.b);
        triangle_buffer.push_back(t.b.x);
        triangle_buffer.push_back(t.b.y);
        triangle_buffer.push_back(t.color.r);
        triangle_buffer.push_back(t.color.g);
        triangle_buffer.push_back(t.color.b);
        triangle_buffer.push_back(t.c.x);
        triangle_buffer.push_back(t.c.y);
        triangle_buffer.push_back(t.color.r);
        triangle_buffer.push_back(t.color.g);
        triangle_buffer.push_back(t.color.b);
    }
    GLuint triangles_vbo = uploadVertices(triangle_buffer.data(), triangle_buffer.size() * sizeof(float));
    // Attrib pointers
    GLuint triangle_vao = createVertexArrayObject();
    glBindVertexArray(triangle_vao);
    setFormat("position 2f color 3f", shaderProgram);

	glm::mat4 proj, view, model;
	int width, height;
	double ratio, timer;

	double mouseX, mouseY;
	
	while (!glfwWindowShouldClose(window))
	{
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);
		// Reupload vertices
		glBindBuffer(GL_ARRAY_BUFFER, polygonBufferRef);
		BufferWriter<float> buffer(polygonBuffer.size());
		b1.addToBuffer(buffer);
		// b2.addToBuffer(buffer);
		for (int i = 0; i < edges.size(); i ++)
		{
            addToBuffer(edges[i], b1, buffer);
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
		//
		glfwGetCursorPos(window, &mouseX, &mouseY);

		// b1 approaches mouse pointer
		float xTarget = mouseX - width/2;
		float yTarget = - mouseY + height/2;
		// b1.velocity().x += (xTarget - b1.position().x) / 100;
		// b1.velocity().y += (yTarget - b1.position().y) / 100;
		// b1.velocity() *= 0.9;

		world.timestep(30);
        b2.shape().color = glm::vec3(0, 0, 0);
		// Check if 'almost' colliding
		if (bodies.overlaps(b1, b2)) {
			// b1.shape().color.r = 1; b1.shape().color.g = 0; b1.shape().color.b = 0;
            b1.shape().color = glm::vec3(1, 0, 0);
		} else {
			// b1.shape().color.r = 1; b1.shape().color.g = 1; b1.shape().color.b = 1;
            b1.shape().color = glm::vec3(0.6f, 0.6f, 0.6f);
		}


		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		ratio = width / (float) height;
		timer = (float)glfwGetTime() * 2;

		glUseProgram(shaderProgram);
			proj = ortho2D(width, height, 0, 1);
			view = viewMatrix2D(0, 0, 1, 1);
			model = glm::mat4{};
			glUniformMatrix4fv(uni_proj, 1, GL_FALSE, glm::value_ptr(proj));
			glUniformMatrix4fv(uni_view, 1, GL_FALSE, glm::value_ptr(view));
			glUniformMatrix4fv(uni_model, 1, GL_FALSE, glm::value_ptr(model));

            glBindVertexArray(triangle_vao);
            glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
            glDrawArrays(GL_TRIANGLES, 0, triangle_buffer.size());
			glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, polygonBufferRef);
            glDrawArrays(GL_LINES, 0, polygonBuffer.size()/2);

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
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	switch(key)
	{
		case GLFW_KEY_UP:
			b1.position().y += 10;
			break;
		case GLFW_KEY_DOWN:
			b1.position().y -= 10;
			break;
		case GLFW_KEY_LEFT:
			b1.position().x -= 10;
			break;
		case GLFW_KEY_RIGHT:
			b1.position().x += 10;
			break;
	}
}
