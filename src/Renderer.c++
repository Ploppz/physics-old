#include <vector>
#include <algorithm>

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

const int lines_vbo_size = 300;

void Renderer::uploadVertices()
{
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    start_indices.clear();
    const int VERT_SIZE = 2; // number of floats

    // Estimate size
    int size = 0;
    for (int i = 0; i < system.numBodies(); i ++)
    {
        size += system.getBody(i).shape().vertices.size() * 2;
    }
    
    BufferWriter<float> buffer(size);
    for (int i = 0; i < system.numBodies(); i ++)
    {
        start_indices.push_back(buffer.getCurrentSize() / VERT_SIZE);
        system.getBody(i).shape().appendStencilTriangles(buffer);
    }
    start_indices.push_back(buffer.getCurrentSize() / VERT_SIZE); // This line is just so we don't need to check end of vector
}

Renderer::Renderer(BodySystem& system)
    : system(system), color1(1, 1, 1), color2(0, 0, 0)
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



    glUseProgram(pos2_program);
    // Triangles VBO
	glGenBuffers(1,&triangles_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
	glBufferData(GL_ARRAY_BUFFER, 3000, NULL, GL_STATIC_DRAW);
    // Triangles VAO
    glGenVertexArrays(1, &triangles_vao);
    glBindVertexArray(triangles_vao);
    setFormat("position 2f", pos2_program);
    // Lines VBO
	glGenBuffers(1,&lines_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
	glBufferData(GL_ARRAY_BUFFER, lines_vbo_size * sizeof(float), NULL, GL_STREAM_DRAW);
    // Lines VAO
    glGenVertexArrays(1, &lines_vao);
    glBindVertexArray(lines_vao);
    setFormat("position 2f", pos2_program);

    glUseProgram(color_program);
    // Red quad
    //color_vbo = uploadVertices(quad, sizeof(quad));
	glGenBuffers(1, &color_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

    glGenVertexArrays(1, &color_vao);
    glBindVertexArray(color_vao);
    setFormat("position 2f", color_program);
}

void Renderer::render(int width, int height, float zoom)
{
    uploadVertices();
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
    assert(start_indices.size() > 0); // The loop still loops when size == 0....... 
    for (int i = 0; i < start_indices.size() - 1; i ++) 
    {
        position = system.getBody(i).real_position();
        length = start_indices[i + 1] - start_indices[i];
        model = glm::translate(glm::mat4{}, glm::vec3(position, 0));
        glUniformMatrix4fv(uni_model, 1, GL_FALSE, glm::value_ptr(model));
        glDrawArrays(GL_TRIANGLE_FAN, start_indices[i], length);
    }

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);        // Enable color buffer
// Draw quad - color polygons
    glUseProgram(color_program);
    glBindVertexArray(color_vao);
    glBindBuffer(GL_ARRAY_BUFFER, color_vbo);

    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);                  // Disable writing to stencil buffer

    glUniform3f(uni_color, color1.r, color1.g, color1.b);
    glStencilFunc(GL_NOTEQUAL, 0, 0xFF);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    glUniform3f(uni_color, color2.r, color2.g, color2.b);
    glStencilFunc(GL_EQUAL, 0, 0xFF);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    glStencilFunc(GL_ALWAYS, 0, 0xFF); // Disable stencil testing
    glDisable(GL_STENCIL_TEST);
// Reupload and draw lines
    glUseProgram(pos2_program);
    glBindVertexArray(lines_vao);
    glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
    model = glm::mat4{};
    glUniformMatrix4fv(uni_model, 1, GL_FALSE, glm::value_ptr(model));
    BufferWriter<float> buffer(lines_vbo_size);
    for (float f : lines_buffer) buffer.write(f);
    // Draw
    glDrawArrays(GL_LINES, 0, lines_buffer.size() / 2);
    lines_buffer.clear(); 
}

void Renderer::setColor1(float r, float g, float b)
{
    color1 = glm::vec3(r, g, b);
}
void Renderer::setColor2(float r, float g, float b)
{
    color2 = glm::vec3(r, g, b);
}

void Renderer::addDot(glm::vec2 dot)
{
    const int radius = 2;

    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y - radius);
    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y + radius);

    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y - radius);
    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y + radius);

    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y + radius);
    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y - radius);

    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y + radius);
    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y - radius);
}
