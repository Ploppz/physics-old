#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>


#include "geometry/Polygon.h"

/* External */
class FontRenderer;
class BodySystem;
class Body;

enum RenderFlags {
    POLYGON_SHOW_VERTEX_NUMBERS = 1,
    POLYGON_SHOW_VELOCITY = 1 <<1,
};

// For now -- decompose once, upload once
class Renderer
{
public:
    Renderer(BodySystem& system);
    void render(float center_x, float center_y, int width, int height, float zoom);
    void set_color_1(float r, float g, float b);
    void set_color_2(float r, float g, float b);
    void set_render_flags(int composite_flag);

    glm::vec2 center_screen_position(float center_x, float center_y, int width, int height, float zoom);

    // Functions to help visualize things.
    void add_dot(glm::vec2 dot);
    void add_vector(glm::vec2 point, glm::vec2 vec);
    void add_polygon_lines(Polygon& p, glm::vec2 position);
    std::vector<float> lines_buffer;
    // Flags
    void set_render_flag(int flag);
    bool has_flag(int flag);


    /** Custom rendering **/
    void append_velocity_lines_to_buffer(Body body);
    void append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer);
    void write_vertex_numbers(Polygon& p);
    void write_distances_to(Polygon& subject, Polygon &other);
    void append_lines_to_vector(Polygon& p);
    void append_lines_to_vector(Intersection& intersection);
private:
    void upload_vertices();
    
    /*** Values ***/
    FontRenderer* font_renderer;
    BodySystem& system;
    std::vector<int> start_indices;
    int render_flags;

    glm::vec3 color1, color2;
    glm::vec3 lines_color;

    /** OpenGL state **/
    // Programs
    GLuint pos2_program, color_program;
    // Uniforms
    GLuint uni_proj, uni_view, uni_center, uni_orientation;
    GLuint uni_color;
    // Uniform block
    /* GLuint model_trans_block_index; */
    // VBO
    GLuint triangles_vbo, color_vbo, lines_vbo;
    // VAO
    GLuint triangles_vao, lines_vao, color_vao;
    /** **/

    // Constants
    const int lines_vbo_size = 6000;
};
