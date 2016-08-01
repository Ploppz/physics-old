#pragma once
#include <vector>
#include <glm/glm.hpp>

#include "glutils.h"
class Polygon;
class Body;
struct Intersection;

class LineBuffer
{
 public:
    void clear_buffer();
    std::vector<float>& get_buffer();
    void write_to_buffer(BufferWriter<float>&);
    void set_color(glm::vec3 color) { this->color = color;}
     /** Render **/
    void write_vertex_numbers(Polygon& p);
    void write_distances_to(Polygon& subject, Polygon &other);
    void append_lines_to_vector(Polygon& p);
    // untransformed
    void append_model_lines_to_vector(Polygon& p);
    void append_lines_to_vector(Intersection& intersection);
    void append_velocity_lines_to_buffer(Body body);
    // Functions to help visualize things.
    void add_dot(glm::vec2 dot);
    void add_dot(glm::vec2 dot, float radius);
    void add_vector(glm::vec2 point, glm::vec2 vec);
    // debugging..
    /* void render_contact(TimeContact c); */
 private:
    bool ok_to_add_floats(int num);
    void add_point(glm::vec2 point);
    void add_point(float x, float y);
 private:
    std::vector<float> buffer;
    glm::vec3 color = glm::vec3(1);
};

