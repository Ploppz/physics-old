#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "Renderer.h"
#include "LineBuffer.h"
#include "../Body.h"
#include "../glutils.h"
#include "../typewriter/FontRenderer.h"
#include "../geometry/Intersection.h"

extern Renderer* g_renderer;

std::vector<float>& LineBuffer::get_buffer() { return buffer; }
void LineBuffer::clear_buffer()
{
    buffer.clear();
}
void LineBuffer::write_to_buffer(BufferWriter<float>& buffer_writer)
{
    for (float f : buffer)
        buffer_writer.write(f);
}



void LineBuffer::write_vertex_numbers(Polygon& p)
{
    for (uint i = 0; i < p.vertices.size(); i ++)
    {
        glm::vec2 point = p.transformed(i);
        g_renderer->get_font_renderer()->setColor(255, 255, 255);
        g_renderer->get_font_renderer()->addText(std::to_string(i), point.x, point.y,  false);
    }
}
void LineBuffer::write_distances_to(Polygon& subject, Polygon &other)
{
    for (uint i = 0; i < subject.vertices.size(); i ++)
    {
        glm::vec2 transformed = subject.transform(subject.vertices[i]);
        int closest_edge;
        float closest_edge_alpha;
        float distance_from_other = distance(transformed, other, closest_edge, closest_edge_alpha);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << distance_from_other;
        g_renderer->get_font_renderer()->setColor(255, 255, 255);
        g_renderer->get_font_renderer()->addText(stream.str(), transformed.x, transformed.y,  false);
    }
}
void LineBuffer::append_lines_to_vector(Polygon& p)
{
    for (uint i = 0; i < p.vertices.size(); i ++) {
        int j = i + 1; j %= p.vertices.size();
        glm::vec2 vec_i = p.transformed(i);
        glm::vec2 vec_j = p.transformed(j);
        buffer.push_back(vec_i.x);
        buffer.push_back(vec_i.y);
        buffer.push_back(vec_j.x);
        buffer.push_back(vec_j.y);
    }
}
void LineBuffer::append_model_lines_to_vector(Polygon& p)
{
    for (uint i = 0; i < p.vertices.size(); i ++) {
        int j = i + 1; j %= p.vertices.size();
        glm::vec2 vec_i = p.vertices[i];
        glm::vec2 vec_j = p.vertices[j];
        buffer.push_back(vec_i.x);
        buffer.push_back(vec_i.y);
        buffer.push_back(vec_j.x);
        buffer.push_back(vec_j.y);
    }
}

void LineBuffer::append_lines_to_vector(Intersection& intersection)
{
    // Loop through HybridVertices
    std::vector<HybridVertex>::iterator next;
    for (auto it = intersection.vertices.begin(); it != intersection.vertices.end(); it ++)
    {
		next = it; next ++;
		if (next == intersection.vertices.end())
            next = intersection.vertices.begin();

        
        glm::vec2 vec_i = it->point;
        glm::vec2 vec_j = next->point;
        buffer.push_back(vec_i.x);
        buffer.push_back(vec_i.y);
        buffer.push_back(vec_j.x);
        buffer.push_back(vec_j.y);
    }
}

void LineBuffer::append_velocity_lines_to_buffer(Body body)
{
    float empty_angle = 1;
    float rotation_radius = 10;

    glm::vec2 center = body.position();
    glm::vec2 velocity = body.velocity();
    float rotation = body.rotation();

    /** Rotational velocity **/
    float velocity_angle = angle_of_vector(velocity);
    bool first_iteration = true;;
    int resolution = 30;
    for (int i = 0; i < resolution; i ++)
    {
        // PERFORMANCE calculate once (more code)
        float angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)i / resolution) * (6.28 - empty_angle));
        glm::vec2 point = center + unit_vector(angle) * rotation_radius;
        float next_angle = velocity_angle + sign(rotation) * (0.5f*empty_angle + ((float)(i+1) / resolution) * (6.28 - empty_angle));
        glm::vec2 next_point = center + unit_vector(next_angle) * rotation_radius;
        if (i < resolution - 1) {
            buffer.push_back(point.x);
            buffer.push_back(point.y);
            buffer.push_back(next_point.x);
            buffer.push_back(next_point.y);
        } else { // Draw vector
            add_vector(point, next_point - point);
        }
    }

    /** Linear velocity **/
    add_vector(center, velocity);
}




void LineBuffer::add_dot(glm::vec2 dot)
{
    const int radius = 1;

    buffer.push_back(dot.x - radius);
    buffer.push_back(dot.y - radius);
    buffer.push_back(dot.x + radius);
    buffer.push_back(dot.y + radius);

    buffer.push_back(dot.x + radius);
    buffer.push_back(dot.y - radius);
    buffer.push_back(dot.x - radius);
    buffer.push_back(dot.y + radius);
}
void LineBuffer::add_vector(glm::vec2 point, glm::vec2 vec)
{
    const int radius = 4;
    const float arrow_angle = 2.4f;

    float vec_angle = atan2(vec.y, vec.x);
    glm::vec2 a1 = glm::vec2(cos(vec_angle - arrow_angle) * radius, sin(vec_angle - arrow_angle) * radius);
    glm::vec2 a2 = glm::vec2(cos(vec_angle + arrow_angle) * radius, sin(vec_angle + arrow_angle) * radius);

    buffer.push_back(point.x);
    buffer.push_back(point.y);
    buffer.push_back(point.x + vec.x);
    buffer.push_back(point.y + vec.y);

    buffer.push_back(point.x + vec.x);
    buffer.push_back(point.y + vec.y);
    buffer.push_back(point.x + vec.x + a1.x);
    buffer.push_back(point.y + vec.y + a1.y);

    buffer.push_back(point.x + vec.x);
    buffer.push_back(point.y + vec.y);
    buffer.push_back(point.x + vec.x + a2.x);
    buffer.push_back(point.y + vec.y + a2.y);

}


/* void LineBuffer::render_contact(TimeContact contact)
{
    glm::vec2 p = contact.ref_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
    p = contact.subj_point.point_t();
    add_dot(p);
    add_vector(p, glm::vec2(contact.normal.x * 10, contact.normal.y * 10));
} */
