/** TODO
- Make depth constants.
**/
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "shaders.h"
#include "tmp.h"
#include "Renderer.h"
#include "BodySystem.h"
#include "geometry/geometry.h"
#include "typewriter/FontTexture.h"
#include "typewriter/FontRenderer.h"
/* glm */
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

GLfloat quad[] = {
	// x, y, tex_x, tex_y
	-1.0f, -1.0f,
	1.0f, -1.0f,	
    1.0f, 1.0f,

	-1.0f, -1.0f,
	1.0f, 1.0f,
	-1.0f, 1.0f
};


void Renderer::upload_vertices()
{
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    start_indices.clear();
    const int VERT_SIZE = 2; // number of floats

    // Estimate size
    int size = 0;
    for (int i = 0; i < system.num_bodies(); i ++)
    {
        size += system.get_body(i).shape().vertices.size() * 2;
    }
    
    BufferWriter<float> buffer(size);
    Body it = system.get_body(0);
    while (it.is_valid())
    {
        start_indices.push_back(buffer.get_current_size() / VERT_SIZE);
        append_stencil_triangle_fan(it.shape(), buffer);
        if (has_flag(POLYGON_SHOW_VERTEX_NUMBERS)) {
            write_vertex_numbers(it.shape());
        }
        if (has_flag(POLYGON_SHOW_VELOCITY)) {
            append_velocity_lines_to_buffer(it); 
        }
        ++ it;
    }
    buffer.unmap();
    start_indices.push_back(buffer.get_current_size() / VERT_SIZE); // This line is just so we don't need to check end of vector
}

Renderer::Renderer(BodySystem& system)
    : system(system), color1(1, 1, 1), color2(0, 0, 0), render_flags {}
{

	/** Typewriter **/
    font_renderer = new FontRenderer(1, "fonts/peep-07x14.bdf");
    font_renderer->setup();

    /** OpenGL **/
    // Create shader program
	GLuint vertexShader, fragmentShader; // unused
    pos2_program  = createShaderProgram(shaders::shaders_pos2_v, shaders::shaders_pos2_f, vertexShader, fragmentShader);
    color_program = createShaderProgram(shaders::shaders_color_v, shaders::shaders_color_f, vertexShader, fragmentShader);
    // Get uniform locations
	uni_proj = glGetUniformLocation(pos2_program, "proj");
	uni_view = glGetUniformLocation(pos2_program, "view");
	uni_center = glGetUniformLocation(pos2_program, "center");
	uni_orientation = glGetUniformLocation(pos2_program, "orientation");
    uni_color = glGetUniformLocation(color_program, "color");
    /* model_trans_block_index = glGetUniformBlockIndex(pos2_program, "ModelTransform"); */



    glUseProgram(pos2_program);
    // Triangles VBO
	glGenBuffers(1,&triangles_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
	glBufferData(GL_ARRAY_BUFFER, 3000, NULL, GL_DYNAMIC_DRAW);
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
	glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_DYNAMIC_DRAW);

    glGenVertexArrays(1, &color_vao);
    glBindVertexArray(color_vao);
    setFormat("position 2f", color_program);
}

glm::vec2 Renderer::center_screen_position(float center_x, float center_y, int width, int height, float zoom)
{
    // note: center_x is the world position that corresponds to center of screen
    glm::mat4 proj = ortho2D(width * zoom, height * zoom, 0, 1);
    glm::mat4 view = viewMatrix2D(center_x, center_y, 1, 1);

    // std::cout << glm::to_string(proj * view) << std::endl;
    glm::mat4 inv_PV = glm::inverse(proj * view);
    glm::vec4 center_of_screen(width / 2, height / 2, 0, 1);
    glm::vec4 result = inv_PV *  center_of_screen;
    return glm::vec2(result.x, result.y);
}
void Renderer::render(float center_x, float center_y, int width, int height, float zoom)
{
    if ( system.last_contact.rewind_time != -1) {
        render_contact(system.last_contact);
    }
    upload_vertices();
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
    view = viewMatrix2D(center_x, center_y, 1, 1);
    glUniformMatrix4fv(uni_proj, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(uni_view, 1, GL_FALSE, glm::value_ptr(view));
    glBindVertexArray(triangles_vao);
    glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    int length;
    glm::vec2 position;
    assert(start_indices.size() > 0); // The loop still loops when size == 0....... 
    for (int i = 0; i < start_indices.size() - 1; i ++) 
    {
        position = system.get_body(i).position();
        length = start_indices[i + 1] - start_indices[i];

        model = glm::translate(glm::mat4{}, glm::vec3(position, 0));
        glUniform2f(uni_center, position.x, position.y);
        glUniform1f(uni_orientation, system.get_body(i).orientation());
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

    glUniform2f(uni_center, 0, 0);
    glUniform1f(uni_orientation, 0);

    //
    const int size_float_buffer = lines_buffer.size() + system.auxilliary_lines.size();
    BufferWriter<float> buffer(size_float_buffer); 
    for (float f : lines_buffer) buffer.write(f);
    for (float f : system.auxilliary_lines) buffer.write(f);
    // Draw
    glDrawArrays(GL_LINES, 0, size_float_buffer / 2);
    lines_buffer.clear(); 
    buffer.unmap();
/** Draw text **/
    font_renderer->render(center_x, center_y, width, height, zoom);
    font_renderer->clearBuffer();
}

void Renderer::set_color_1(float r, float g, float b)
{
    color1 = glm::vec3(r, g, b);
}
void Renderer::set_color_2(float r, float g, float b)
{
    color2 = glm::vec3(r, g, b);
}
bool Renderer::has_flag(int flag)
{
    return (render_flags & flag) != 0;
}
void Renderer::set_render_flag(int flag)
{
    render_flags |= flag;
}

void Renderer::add_dot(glm::vec2 dot)
{
    const int radius = 1;

    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y - radius);
    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y + radius);

    lines_buffer.push_back(dot.x + radius);
    lines_buffer.push_back(dot.y - radius);
    lines_buffer.push_back(dot.x - radius);
    lines_buffer.push_back(dot.y + radius);
}
void Renderer::add_vector(glm::vec2 point, glm::vec2 vec)
{
    const int radius = 4;
    const float arrow_angle = 2.4f;

    float vec_angle = atan2(vec.y, vec.x);
    glm::vec2 a1 = glm::vec2(cos(vec_angle - arrow_angle) * radius, sin(vec_angle - arrow_angle) * radius);
    glm::vec2 a2 = glm::vec2(cos(vec_angle + arrow_angle) * radius, sin(vec_angle + arrow_angle) * radius);

    lines_buffer.push_back(point.x);
    lines_buffer.push_back(point.y);
    lines_buffer.push_back(point.x + vec.x);
    lines_buffer.push_back(point.y + vec.y);

    lines_buffer.push_back(point.x + vec.x);
    lines_buffer.push_back(point.y + vec.y);
    lines_buffer.push_back(point.x + vec.x + a1.x);
    lines_buffer.push_back(point.y + vec.y + a1.y);

    lines_buffer.push_back(point.x + vec.x);
    lines_buffer.push_back(point.y + vec.y);
    lines_buffer.push_back(point.x + vec.x + a2.x);
    lines_buffer.push_back(point.y + vec.y + a2.y);

}


void Renderer::append_velocity_lines_to_buffer(Body body)
{
    float empty_angle = 1;
    float rotation_radius = 10;

    glm::vec2 center = body.position();
    glm::vec2 velocity = body.velocity();
    float rotation = body.rotation();

    /** Rotational velocity **/
    float velocity_angle = angle_of_vector(velocity);
    bool first_iteration = true;;
    int resolution = 30;
    for (int i = 0; i < resolution; i ++)
    {
        // PERFORMANCE calculate once (more code)
        float angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)i / resolution) * (6.28 - empty_angle));
        glm::vec2 point = center + unit_vector(angle) * rotation_radius;
        float next_angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)(i+1) / resolution) * (6.28 - empty_angle));
        glm::vec2 next_point = center + unit_vector(next_angle) * rotation_radius;
        if (i < resolution - 1) {
            lines_buffer.push_back(point.x);
            lines_buffer.push_back(point.y);
            lines_buffer.push_back(next_point.x);
            lines_buffer.push_back(next_point.y);
        } else { // Draw vector
            add_vector(point, next_point - point);
        }
    }

    /** Linear velocity **/
    add_vector(center, velocity);
}


/** Polygon **/

void Renderer::append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer)
{
    for (uint i = 0; i < p.vertices.size(); i ++)
    {
        buffer.write(p.vertices[i].x, p.vertices[i].y);
    }
}
void Renderer::write_vertex_numbers(Polygon& p)
{
    for (uint i = 0; i < p.vertices.size(); i ++)
    {
        glm::vec2 point = p.transformed(i);
        font_renderer->setColor(255, 255, 255);
        font_renderer->addText(std::to_string(i), point.x, point.y,  false);
    }
}
void Renderer::write_distances_to(Polygon& subject, Polygon &other)
{
    for (uint i = 0; i < subject.vertices.size(); i ++)
    {
        glm::vec2 transformed = subject.transform(subject.vertices[i]);
        int closest_edge;
        float closest_edge_alpha;
        float distance_from_other = distance(transformed, other, closest_edge, closest_edge_alpha);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << distance_from_other;
        font_renderer->setColor(255, 255, 255);
        font_renderer->addText(stream.str(), transformed.x, transformed.y,  false);
    }
}
void Renderer::append_lines_to_vector(Polygon& p)
{
    for (uint i = 0; i < p.vertices.size(); i ++) {
        int j = i + 1; j %= p.vertices.size();
        glm::vec2 vec_i = p.transform(p.vertices[i]);
        glm::vec2 vec_j = p.transform(p.vertices[j]);
        lines_buffer.push_back(vec_i.x);
        lines_buffer.push_back(vec_i.y);
        lines_buffer.push_back(vec_j.x);
        lines_buffer.push_back(vec_j.y);
    }
}


/** Intersection **/


void Renderer::append_lines_to_vector(Intersection& intersection)
{
    // Loop through HybridVertices
    std::vector<HybridVertex>::iterator next;
    for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
    {
		next = it; next ++;
		if (next == intersection.vertices.end())
            next = intersection.vertices.begin();

        
        glm::vec2 vec_i = it->point;
        glm::vec2 vec_j = next->point;
        lines_buffer.push_back(vec_i.x);
        lines_buffer.push_back(vec_i.y);
        lines_buffer.push_back(vec_j.x);
        lines_buffer.push_back(vec_j.y);
    }
}
void Renderer::render_contact(Contact contact)
{
    glm::vec2 p = contact.ref_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
    p = contact.subj_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
}
