#include <vector>

#include "shaders.h"
#include "tmp.h"
#include "Renderer.h"
#include "BodySystem.h"
#include "Geometry.h"
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
    pos2_program  = createShaderProgram(shaders::shaders_pos2_v, shaders::shaders_pos2_f, vertexShader, fragmentShader);
    color_program = createShaderProgram(shaders::shaders_color_v, shaders::shaders_color_f, vertexShader, fragmentShader);
    // Get uniform locations
	uni_proj = glGetUniformLocation(pos2_program, "proj");
	uni_view = glGetUniformLocation(pos2_program, "view");
	uni_model = glGetUniformLocation(pos2_program, "model");
    uni_color = glGetUniformLocation(color_program, "color");

    // Decompose polygons
    std::vector<float> triangle_buffer;
    for (int i = 0; i < system.numBodies(); i ++)
    {
        start_indices.push_back(triangle_buffer.size() / 2);
        system.getBody(i).shape().appendStencilTriangles(triangle_buffer);
    }
    start_indices.push_back(triangle_buffer.size() / 2); // This line is just so we don't need to check end of vector

    glUseProgram(pos2_program);
    // Upload triangles
    // triangles_vbo = uploadVertices(triangle_buffer.data(), triangle_buffer.size() * sizeof(float));
	glGenBuffers(1,&triangles_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
	glBufferData(GL_ARRAY_BUFFER, triangle_buffer.size() * sizeof(float), triangle_buffer.data(), GL_STATIC_DRAW);
    // Attrib pointers
    triangles_vao = createVertexArrayObject();
    glBindVertexArray(triangles_vao);
    setFormat("position 2f", pos2_program);

    glUseProgram(color_program);
    // Red quad
    color_vbo = uploadVertices(quad, sizeof(quad));
    color_vao = createVertexArrayObject();
    glBindVertexArray(color_vao);
    setFormat("position 2f", color_program);
}

void Renderer::render(int width, int height, float zoom)
{
    glClear(GL_STENCIL_BUFFER_BIT);
    glEnable(GL_STENCIL_TEST);
// Draw stencil triangles
    glUseProgram(pos2_program);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glStencilFunc(GL_ALWAYS, 0, 0xFF);
    glStencilMask(0xFF);
    glStencilOp(GL_INVERT, GL_INVERT, GL_INVERT);
	glm::mat4 proj, view, model;
    proj = ortho2D(width * zoom, height * zoom, 0, 1);
    view = viewMatrix2D(0, 0, 1, 1);
    glUniformMatrix4fv(uni_proj, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(uni_view, 1, GL_FALSE, glm::value_ptr(view));
    glBindVertexArray(triangles_vao);
    glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    int length;
    glm::vec2 position;
    for (int i = 0; i < start_indices.size() - 1; i ++) 
    {
        position = system.getBody(i).real_position();
        length = start_indices[i + 1] - start_indices[i];
        model = glm::translate(glm::mat4{}, glm::vec3(position, 0));
        glUniformMatrix4fv(uni_model, 1, GL_FALSE, glm::value_ptr(model));
        glDrawArrays(GL_TRIANGLE_FAN, start_indices[i], length);
    }

// Draw quad - color polygons
    glUseProgram(color_program);
    glBindVertexArray(color_vao);
    glUniform3f(uni_color, 0, 0.6f, 0);
    glBindBuffer(GL_ARRAY_BUFFER, color_vbo);

    glStencilFunc(GL_NOTEQUAL, 0, 0xFF);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}
