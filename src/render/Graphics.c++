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
#include "Graphics.h"
#include "LineRenderer.h"
#include "font/FontTexture.h"
#include "font/FontRenderer.h"

#include "shaders.h"
#include "tmp.h"
#include "BodySystem.h"
#include "Body.h"
#include "geometry/geometry.h"
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



Graphics::Graphics(BodySystem& system)
    : system(system), render_flags {}, color1(1, 1, 1), color2(0, 0, 0)
{
    /** Set up buffers **/
    set_buffers_to_default();

    /** OpenGL **/
    // Create shader program
    pos2_program  = createShaderProgram(shaders::shaders_pos2_v, shaders::shaders_pos2_f);
    color_program = createShaderProgram(shaders::shaders_color_v, shaders::shaders_color_f);
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
	glBufferData(GL_ARRAY_BUFFER, TRIANGLE_FAN_VBO_SIZE, NULL, GL_DYNAMIC_DRAW);
    // Triangles VAO
    glGenVertexArrays(1, &triangles_vao);
    glBindVertexArray(triangles_vao);
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

void Graphics::set_buffers_to_default()
{
    set_active_buffers("default");
}

void Graphics::set_active_buffers(const std::string& name)
{
    DebugBegin();
    if (std::find(buffer_names.begin(), buffer_names.end(), name) == buffer_names.end())
        buffer_names.push_back(name);

    line_renderer.set_active_buffer(name);
    font_renderer.set_active_buffer(name);
}

void Graphics::upload_vertices()
{
	glBindBuffer(GL_ARRAY_BUFFER, triangles_vbo);
    start_indices.clear();
    const int VERT_SIZE = 2; // number of floats

    // Estimate size
    int size = 0;
    for (int i = 0; i < system.num_bodies(); i ++)
    {
        size += system.get_body(i).shape().num_vertices() * 2;
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
            append_velocity_lines(it); 
        }
        ++ it;
    }
    buffer.unmap();
    start_indices.push_back(buffer.get_current_size() / VERT_SIZE); // This line is just so we don't need to check end of vector
}

glm::vec2 Graphics::center_screen_position(float center_x, float center_y, int width, int height, float zoom)
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

void Graphics::render_classic(float center_x, float center_y, int width, int height, float zoom)
{
    upload_vertices();
    glClear(GL_STENCIL_BUFFER_BIT);
    glEnable(GL_STENCIL_TEST);
// Draw stencil triangles //
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
// Draw quad - color polygons //
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


}
void Graphics::render(float center_x, float center_y, int width, int height, float zoom)
{
    line_renderer.render(center_x, center_y, width, height, zoom);
    font_renderer.render(center_x, center_y, width, height, zoom);
}
void Graphics::render_all_buffers(float center_x, float center_y, int width, int height, float zoom)
{
    DebugBegin();
    for (const std::string& buffer_name : buffer_names)
    {
        set_active_buffers(buffer_name);
        render(center_x, center_y, width, height, zoom);

        // Because the default buffers should be renewed every frame:
        if (buffer_name == "default")
            clear_buffers();
    }
}
void Graphics::render(StatisticsCollection statistics)
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

void Graphics::clear_buffers()
{
    // TODO maybe remove (from public at least) since we should do this pretty much only when rendering
    line_renderer.clear_buffer();
    font_renderer.clear_buffer();
}

void Graphics::set_color_1(float r, float g, float b) {
    color1 = glm::vec3(r, g, b);
}

void Graphics::set_color_2(float r, float g, float b) {
    color2 = glm::vec3(r, g, b);
}

void Graphics::set_line_color(float r, float g, float b) {
    line_renderer.set_color(r, g, b);
}

void Graphics::set_font_color(float r, float g, float b) {
    font_renderer.set_color(r, g, b);
}

bool Graphics::has_flag(int flag) {
    return (render_flags & flag) != 0;
}

void Graphics::set_render_flag(int flag) {
    render_flags |= flag;
}






/** Polygon **/

// TODO: just an example of what we should move to a StencilRenderer
void Graphics::append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer)
{
    for (Vertex v : p.model_vertices())
    {
        buffer.write(v.point.x, v.point.y);
    }
}


//////////////////////////////////////////////////////
//
// Draw commands
//
/////////////////////////////////////////////////////


void Graphics::write_vertex_numbers(Polygon& p)
{
    for (Vertex v : p.vertices())
    {
        font_renderer.set_color(1,1,1);
        font_renderer.add_text(std::to_string(v.index), v.point.x, v.point.y,  false);
    }
}
void Graphics::write_distances_to(Polygon& subject, Polygon &other)
{
    for (Vertex v : subject.vertices())
    {
        int closest_edge;
        float closest_edge_alpha;
        float distance_from_other = distance(v.point, other, closest_edge, closest_edge_alpha);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << distance_from_other;
        font_renderer.set_color(255, 255, 255);
        font_renderer.add_text(stream.str(), v.point.x, v.point.y,  false);
    }
}
void Graphics::append_lines(Polygon& p)
{
    for (Edge e : p.edges())
        line_renderer.add_line(e.start, e.end);
}
void Graphics::append_model_lines(Polygon& p)
{
    for (Edge e : p.model_edges())
        line_renderer.add_line(e.start, e.end);
}

void Graphics::append_lines(Intersection& intersection)
{
    // Loop through HybridVertices
    std::vector<HybridVertex>::iterator next;
    for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
    {
		next = it; next ++;
		if (next == intersection.vertices.end())
            next = intersection.vertices.begin();

        line_renderer.add_line(it->point, next->point);
    }
}

void Graphics::append_velocity_lines(Body body)
{
    float empty_angle = 1;
    float rotation_radius = 10;

    glm::vec2 center = body.position();
    glm::vec2 velocity = body.velocity();
    float rotation = body.rotation();

    /** Rotational velocity **/
    float velocity_angle = angle_of_vector(velocity);
    bool first_iteration = true;
    int resolution = 30;
    for (int i = 0; i < resolution; i ++)
    {
        // PERFORMANCE calculate once (more code)
        float angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)i / resolution) * (6.28 - empty_angle));
        glm::vec2 point = center + unit_vector(angle) * rotation_radius;
        float next_angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)(i+1) / resolution) * (6.28 - empty_angle));
        glm::vec2 next_point = center + unit_vector(next_angle) * rotation_radius;
        if (i < resolution - 1) {
            line_renderer.add_line(point, next_point);
        } else { // Draw vector
            line_renderer.add_vector(point, next_point - point);
        }
    }

    /** Linear velocity **/
    line_renderer.add_vector(center, velocity);
}



void Graphics::labelled_line(glm::vec2 start, glm::vec2 end, const std::string& label)
{
    line_renderer.add_line(start, end);
    glm::vec2 point = (start + end) / 2.f;
    font_renderer.add_text(label, point.x, point.y, 4, false);
}




/* void LineRenderer::render_contact(TimeContact contact)
{
    glm::vec2 p = contact.ref_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
    p = contact.subj_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
} */
