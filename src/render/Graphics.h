#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>


/* src */
#include "LineRenderer.h"
#include "font/FontRenderer.h"

#include "geometry/Polygon.h"
#include "debug/StatisticsCollection.h"

/* External */
class BodySystem;
class Body;
struct TimeContact;

enum RenderFlags {
    POLYGON_SHOW_VERTEX_NUMBERS = 1,
    POLYGON_SHOW_VELOCITY = 1 <<1,
};


//// About Graphics ////
//
// Earlier: Renderer
// It now acts primarily as the interface to several Renderers, such as LineRenderer and FontRenderer.
// It still does some of the rendering itself: such as triangles & stencil rendering of polygons.
// This should change.
// Note to self: The Renderer for polygons should take single polygons.
//   - or should it just draw the basic geometry needed? YES.
//
////////////////////////

class Graphics
{
// Methods //
 public:
    Graphics(BodySystem& system);
    void render(float center_x, float center_y, int width, int height, float zoom);
    void render_classic(float center_x, float center_y, int width, int height, float zoom);
    void render_all_buffers(float center_x, float center_y, int width, int height, float zoom);
    void render(StatisticsCollection statistics);
    void clear_buffers(); // currently only clears line_renderer and font_renderer //
    void set_color_1(float r, float g, float b);
    void set_color_2(float r, float g, float b);
    void set_line_color(glm::vec3 col) { set_line_color(col.r, col.g, col.b); }
    void set_line_color(float r, float g, float b);
    void set_font_color(glm::vec3 col) { set_font_color(col.r, col.g, col.b); }
    void set_font_color(float r, float g, float b);
    void set_render_flags(int composite_flag);

    glm::vec2 center_screen_position(float center_x, float center_y, int width, int height, float zoom);

    /** Flags **/
    void set_render_flag(int flag);
    bool has_flag(int flag);

    /** Buffer management **/
    // applies the buffer with name `name` to all buffer types (lines, fonts...)
    void set_active_buffers(const std::string& name);
    void set_buffers_to_default();

    /** Drawing interface **/
    // Compound object //
    void write_vertex_numbers(Polygon& p);
    void write_distances_to(Polygon& subject, Polygon &other);
    void append_lines(Polygon& p);
    // .. untransformed //
    void append_model_lines(Polygon& p);
    void append_lines(Intersection& intersection);
    void append_velocity_lines(Body body);
    // other
    void labelled_line(glm::vec2 start, glm::vec2 end, const std::string& label);
 private:
    void upload_vertices();
    /** Rendering **/
    void append_stencil_triangle_fan(Polygon& p, BufferWriter<float>& buffer);

// Members //
 public:
    // these are for accessing a more basic interface
    LineRenderer line_renderer;
    FontRenderer font_renderer;
 private:
    BodySystem& system; // TODO unnecessary if we make more modules
    std::vector<std::string> buffer_names; // used just for looping through the buffers

    std::vector<int> start_indices; // start indices of the uploaded stencil triangle fans
    int render_flags;

    glm::vec3 color1, color2;
    glm::vec3 lines_color;

    /* Constants */
    const int LINES_VBO_SIZE = 40000;
    const int TRIANGLE_FAN_VBO_SIZE = 10000;

    /* OpenGL state */
    // Programs
    GLuint pos2_program, pos2col3_program, color_program;
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
