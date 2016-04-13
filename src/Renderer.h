#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>


#include "BodySystem.h"
#include "geometry/Polygon.h"

// For now -- decompose once, upload once
class Renderer
{
public:
    Renderer(BodySystem& system);
    void render(float center_x, float center_y, int width, int height, float zoom);
    void set_color_1(float r, float g, float b);
    void set_color_2(float r, float g, float b);

    glm::vec2 center_screen_position(float center_x, float center_y, int width, int height, float zoom);

    // Functions to help visualize things.
    void add_dot(glm::vec2 dot);
    void add_vector(glm::vec2 point, glm::vec2 vec);
    void add_polygon_lines(Polygon& p, glm::vec2 position);
    std::vector<float> lines_buffer;
private:
    void upload_vertices();
    BodySystem& system;
    std::vector<int> start_indices;

    glm::vec3 color1, color2;
    glm::vec3 lines_color;

    // Programs
    GLuint pos2_program, color_program;
    // Uniforms
    GLuint uni_proj, uni_view, uni_model;
    GLuint uni_color;
    // VBO
    GLuint triangles_vbo, color_vbo, lines_vbo;
    // VAO
    GLuint triangles_vao, lines_vao, color_vao;
};
