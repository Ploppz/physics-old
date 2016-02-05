#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>


#include "BodySystem.h"
#include "Polygon.h"

// For now -- decompose once, upload once
class Renderer
{
public:
    Renderer(BodySystem& system);
    void render(int width, int height, float zoom);
    void setColor1(float r, float g, float b);
    void setColor2(float r, float g, float b);

    // Functions to help visualize things.
    void addDot(glm::vec2 dot);
    void addPolygonLines(Polygon& p, glm::vec2 position);
    std::vector<float> lines_buffer;
private:
    void uploadVertices();
    BodySystem& system;
    std::vector<int> start_indices;

    glm::vec3 color1, color2;

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
