#include <vector>

#include "shaders.h"
#include "tmp.h"
#include "Renderer.h"
#include "BodySystem.h"
/* glm */
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

GLfloat quad[] = {
	// x, y, tex_x, tex_y
	-1.0f, -1.0f,
	1.0f, -1.0f,	
    1.0f, 1.0f,

	-1.0f, -1.0f,
	1.0f, 1.0f,
	-1.0f, 1.0f
};

Renderer::Renderer(BodySystem& system)
    : system(system)
{
    // Create shader program
	GLuint vertexShader, fragmentShader; // unused
    pos2col3_program  = createShaderProgram(shaders::shaders_pos2col3_v, shaders::shaders_pos2col3_f, vertexShader, fragmentShader);
    red_program = createShaderProgram(shaders::shaders_red_v, shaders::shaders_red_f, vertexShader, fragmentShader);
    // Get uniform locations
	uni_proj = glGetUniformLocation(pos2col3_program, "proj");
	uni_view = glGetUniformLocation(pos2col3_program, "view");
	uni_model = glGetUniformLocation(pos2col3_program, "model");

    // Decompose polygons
    std::vector<Triangle> triangles;
	std::vector<LineSegment> lines;
    for (int i = 0; i < system.numBodies(); i ++)
    {
        start_indices.push_back(triangles.size() * 3);
        std::cout << triangles.size() << " TRIANGLES SIZE "<< std::endl;
        system.getBody(i).shape().decompose(triangles, lines);   
    }
    start_indices.push_back(triangles.size() * 3); // This line is just so we don't need to check end of vector

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
    std::cout << "NUM INDICES " << triangle_buffer.size() / 3 << std::endl;


    glUseProgram(pos2col3_program);
    // Upload triangles
    triangles_vbo = uploadVertices(triangle_buffer.data(), triangle_buffer.size() * sizeof(float));
    // Attrib pointers
    triangles_vao = createVertexArrayObject();
    glBindVertexArray(triangles_vao);
    setFormat("position 2f color 3f", pos2col3_program);

    glUseProgram(red_program);
    // Red quad
    red_vbo = uploadVertices(quad, sizeof(quad));
    red_vao = createVertexArrayObject();
    glBindVertexArray(red_vao);
    setFormat("position 2f", red_program);
}

void Renderer::render(int width, int height)
{
    glClear(GL_STENCIL_BUFFER_BIT);
    glEnable(GL_STENCIL_TEST);

    glUseProgram(pos2col3_program);
    glStencilFunc(GL_ALWAYS, 0, 0xFF);
    glStencilMask(0xFF);
    glStencilOp(GL_INCR, GL_INCR, GL_INCR);
	glm::mat4 proj, view, model;
    proj = ortho2D(width, height, 0, 1);
    view = viewMatrix2D(0, 0, 1, 1);
    glUniformMatrix4fv(uni_proj, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(uni_view, 1, GL_FALSE, glm::value_ptr(view));
    glBindVertexArray(triangles_vao);
    glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    int length;
    std::cout << "render" << std::endl;
    for (int i = 0; i < start_indices.size() - 1; i ++) 
    {
        length = start_indices[i + 1] - start_indices[i];
        std::cout << "START LENGTH " << start_indices[i] << " -- " << length << std::endl;
        model = glm::translate(glm::mat4{}, glm::vec3(system.getBody(i).position(), 0));
        glUniformMatrix4fv(uni_model, 1, GL_FALSE, glm::value_ptr(model));
        glDrawArrays(GL_TRIANGLES, start_indices[i], length);
    }

    // Draw quad - color where intersection
    glUseProgram(red_program);
    glBindVertexArray(red_vao);
    glBindBuffer(GL_ARRAY_BUFFER, red_vbo);
    glStencilFunc(GL_EQUAL, 2, 0xFF);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}
