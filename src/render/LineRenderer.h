#pragma once
#include <vector>
#include <glm/glm.hpp>

#include "glutils.h"
#include "Renderer.h"

class Polygon;
class Body;
struct Intersection;

class LineRenderer : public Renderer
{
 public:
    LineRenderer();
    void render(float center_x, float center_y, int width, int height, float zoom);
    void set_color(glm::vec3 color) { this->color = color;}
    void set_color(float r, float g, float b) { set_color(glm::vec3(r,g,b)); }
    // Draw //
    void add_dot(glm::vec2 dot);
    void add_dot(glm::vec2 dot, float radius);
    void add_line(glm::vec2 start, glm::vec2 end);
    void add_vector(glm::vec2 point, glm::vec2 vec);
    void add_aabb(float min_x, float max_x, float min_y, float max_y);
 private:
    bool ok_to_add_floats(int num);
    void add_point(glm::vec2 point);
    void add_point(float x, float y);
 private:
    glm::vec3 color = glm::vec3(1);

    // GL //
    GLuint pos2col3_prg;
    GLuint vbo;
    GLuint vao;
    GLuint uni_proj, uni_view, uni_center, uni_orientation;
    // - constants
    const int VBO_SIZE = 40000; // number of floats
};

