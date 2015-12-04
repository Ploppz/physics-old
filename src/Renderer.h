#pragma once
#include <GL/glew.h>
#include <vector>


#include "BodySystem.h"

// For now -- decompose once, upload once
class Renderer
{
public:
    Renderer(BodySystem& system);
    void render(int width, int height, float zoom);
private:
    BodySystem& system;
    std::vector<int> start_indices;

    // Programs
    GLuint pos2_program, color_program;
    // Uniforms
    GLuint uni_proj, uni_view, uni_model;
    GLuint uni_color;
    // VBO
    GLuint triangles_vbo, color_vbo;
    // VAO
    GLuint triangles_vao, color_vao;
};
