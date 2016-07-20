#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <glm/glm.hpp>

#include "Renderer.h"
#include "LineRenderer.h"
#include "Body.h"
#include "glutils.h"
#include "geometry/Intersection.h"
#include "shaders.h"


const int GL_VERTEX_SIZE = 5;

LineRenderer::LineRenderer()
{
    pos2col3_prg  = createShaderProgram(shaders::shaders_pos2col3_v, shaders::shaders_pos2col3_f);
	uni_proj = glGetUniformLocation(pos2col3_prg, "proj");
	uni_view = glGetUniformLocation(pos2col3_prg, "view");
	uni_center = glGetUniformLocation(pos2col3_prg, "center");
	uni_orientation = glGetUniformLocation(pos2col3_prg, "orientation");
    glUseProgram(pos2col3_prg);
// VBO //
	glGenBuffers(1,&vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, VBO_SIZE * sizeof(float), NULL, GL_STREAM_DRAW);
// VAO //
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    setFormat("position 2f color 3f", pos2col3_prg);
}

void LineRenderer::render(float center_x, float center_y, int width, int height, float zoom)
{
    glUseProgram(pos2col3_prg);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
// Uniforms //
    glUniform2f(uni_center, 0, 0);
    glUniform1f(uni_orientation, 0);
    glm::mat4 proj, view, model;
    proj = ortho2D(width * zoom, height * zoom, 0, 1);
    view = viewMatrix2D(center_x, center_y, 1, 1);
    glUniformMatrix4fv(uni_proj, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(uni_view, 1, GL_FALSE, glm::value_ptr(view));
// Reupload //
    int size_float_buffer = 0;
    for (auto it : buffers)
        size_float_buffer += it.second.size();
    size_float_buffer = std::min(size_float_buffer, VBO_SIZE);

    if (size_float_buffer > 0) {
        BufferWriter<float> buffer(size_float_buffer); 
        for (auto it : buffers)
            buffer.write(it.second);
        buffer.unmap();
    }

// Draw //
    glDrawArrays(GL_LINES, 0, size_float_buffer / GL_VERTEX_SIZE); // 5: maybe make it a constant?

}



void LineRenderer::add_point(glm::vec2 point)
{
    active_buffer->push_back(point.x);
    active_buffer->push_back(point.y);
    active_buffer->push_back(color.r);
    active_buffer->push_back(color.g);
    active_buffer->push_back(color.b);
}

void LineRenderer::add_point(float x, float y)
{
    active_buffer->push_back(x);
    active_buffer->push_back(y);
    active_buffer->push_back(color.r);
    active_buffer->push_back(color.g);
    active_buffer->push_back(color.b);
}


void LineRenderer::add_dot(glm::vec2 dot, float radius) 
{
    add_point(dot - glm::vec2(radius));
    add_point(dot + glm::vec2(radius));
    add_point(dot.x + radius, dot.y - radius);
    add_point(dot.x - radius, dot.y + radius);
}

void LineRenderer::add_dot(glm::vec2 dot)
{
    add_dot(dot, 2);
}

void LineRenderer::add_line(glm::vec2 start, glm::vec2 end)
{
    add_point(start);
    add_point(end);
}

void LineRenderer::add_vector(glm::vec2 point, glm::vec2 vec)
{
    // const float radius = 0.07f;
    const float radius = glm::length(vec) / 15.f;
    const float arrow_angle = 2.7f;

    float vec_angle = atan2(vec.y, vec.x);
    glm::vec2 a1 = glm::vec2(cos(vec_angle - arrow_angle) * radius, sin(vec_angle - arrow_angle) * radius);
    glm::vec2 a2 = glm::vec2(cos(vec_angle + arrow_angle) * radius, sin(vec_angle + arrow_angle) * radius);

    add_point(point);
    add_point(point + vec);;

    add_point(point + vec);
    add_point(point + vec + a1);

    add_point(point + vec);
    add_point(point + vec + a2);
}

void LineRenderer::add_aabb(float min_x, float max_x, float min_y, float max_y)
{
    add_point(min_x, min_y);
    add_point(min_x, max_y);

    add_point(min_x, max_y);
    add_point(max_x, max_y);

    add_point(max_x, max_y);
    add_point(max_x, min_y);

    add_point(max_x, min_y);
    add_point(min_x, min_y);
}
