#pragma once
#include <GL/glew.h>
#include <vector>


#include "BodySystem.h"

// For now -- decompose once, upload once
class Renderer
{
public:
    Renderer(BodySystem& system);
    void render(int, int);
private:
    BodySystem& system;
    std::vector<int> start_indices;

    // Programs
    GLuint pos2col3_program, red_program;
    // Uniforms
    GLuint uni_proj, uni_view, uni_model;
    // VBO
    GLuint triangles_vbo, red_vbo;
    // VAO
    GLuint triangles_vao, red_vao;
};
