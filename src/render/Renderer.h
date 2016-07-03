#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>


/* src */
#include "LineBuffer.h"
#include "geometry/Polygon.h"
#include "debug/StatisticsCollection.h"

/* External */
class FontRenderer;
class BodySystem;
class Body;
struct TimeContact;

enum RenderFlags {
    POLYGON_SHOW_VERTEX_NUMBERS = 1,
    POLYGON_SHOW_VELOCITY = 1 <<1,
};

// For now -- decompose once, upload once
class Renderer
{
/** METHODS **/
 public:
    Renderer(BodySystem& system);
    void render(float center_x, float center_y, int width, int height, float zoom);
    void render(StatisticsCollection statistics);
    void set_color_1(float r, float g, float b);
    void set_color_2(float r, float g, float b);
    void set_render_flags(int composite_flag);
    FontRenderer* get_font_renderer();

    glm::vec2 center_screen_position(float center_x, float center_y, int width, int height, float zoom);

    /** Flags **/
    void set_render_flag(int flag);
    bool has_flag(int flag);


 private:
    void upload_vertices();
    /** Rendering **/
    void append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer);

/** MEMBERS **/
 public:
    LineBuffer extra_line_buffer;
 private:
    BodySystem& system;
    FontRenderer* font_renderer;
    LineBuffer line_buffer;
    std::vector<int> start_indices; // start indices of the uploaded stencil triangle fans
    int render_flags;

    glm::vec3 color1, color2;
    glm::vec3 lines_color;

    /* Constants */
    const int LINES_VBO_SIZE = 40000;

    /* OpenGL state */
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
};
