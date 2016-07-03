/** TODO
- Make depth constants.
**/

#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

/* src */
#include "Renderer.h"
#include "shaders.h"
#include "tmp.h"
#include "BodySystem.h"
#include "Body.h"
#include "geometry/geometry.h"
#include "typewriter/FontTexture.h"
#include "typewriter/FontRenderer.h"
#include "imgui/imgui.h"
#include "debug/StatisticsCollection.h"
#include "debug/debug.h"

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
            line_buffer.write_vertex_numbers(it.shape());
        }
        if (has_flag(POLYGON_SHOW_VELOCITY)) {
            line_buffer.append_velocity_lines_to_buffer(it); 
        }
        ++ it;
    }
    buffer.unmap();
    start_indices.push_back(buffer.get_current_size() / VERT_SIZE); // This line is just so we don't need to check end of vector
}

Renderer::Renderer(BodySystem& system)
    : system(system), render_flags {}, color1(1, 1, 1), color2(0, 0, 0)
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
	glBufferData(GL_ARRAY_BUFFER, LINES_VBO_SIZE * sizeof(float), NULL, GL_STREAM_DRAW);
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
    line_buffer.clear_buffer();
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

    for (int i = 0; i < static_cast<int>(start_indices.size()) - 1; i ++) 
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

    const int size_float_buffer = line_buffer.get_buffer().size() + extra_line_buffer.get_buffer().size();
    if (size_float_buffer > 0) {
        BufferWriter<float> buffer(size_float_buffer); 
        line_buffer.write_to_buffer(buffer);
        extra_line_buffer.write_to_buffer(buffer);
        buffer.unmap();
    }


    // Draw
    glDrawArrays(GL_LINES, 0, size_float_buffer / 2);
/** Draw text **/
    font_renderer->render(center_x, center_y, width, height, zoom);
    font_renderer->clearBuffer();
}
void Renderer::render(StatisticsCollection statistics)
{
    DebugBegin();
    ImGui::CollapsingHeader("Statistics");
    ImGui::Text("\n");
        ImGui::Columns(4, "mycolumns");
        ImGui::Separator();
        ImGui::Text("name"); ImGui::NextColumn();
        ImGui::Text("avg"); ImGui::NextColumn();
        ImGui::Text("min"); ImGui::NextColumn();
        ImGui::Text("max"); ImGui::NextColumn();
        ImGui::Separator();
        for (auto& it : statistics)
        {
            ImGui::Text("%s", it.first.c_str()); ImGui::NextColumn();
            ImGui::Text("%f", it.second.get_avg()); ImGui::NextColumn();
            ImGui::Text("%f", it.second.get_min()); ImGui::NextColumn();
            ImGui::Text("%f", it.second.get_max()); ImGui::NextColumn();
            ImGui::Separator();
        }
        ImGui::Columns(1);

    ImGui::Text("\n");
        ImGui::Columns(2, "mycolumns");
        ImGui::Separator();
        ImGui::Text("name"); ImGui::NextColumn();
        ImGui::Text("count"); ImGui::NextColumn();
        ImGui::Separator();
        for (auto it = statistics.begin_counters(); it != statistics.end_counters(); ++ it)
        {
            ImGui::Text("%s", it->first.c_str()); ImGui::NextColumn();
            ImGui::Text("%d", it->second); ImGui::NextColumn();
            ImGui::Separator();
        }
        ImGui::Columns(1);
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
FontRenderer* Renderer::get_font_renderer()
{
    return font_renderer;
}





/** Polygon **/

void Renderer::append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer)
{
    for (uint i = 0; i < p.vertices.size(); i ++)
    {
        buffer.write(p.vertices[i].x, p.vertices[i].y);
    }
}
