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
    for (Vertex v : p.vertices())
    {
        g_renderer->get_font_renderer()->setColor(255, 255, 255);
        g_renderer->get_font_renderer()->addText(std::to_string(v.index), v.point.x, v.point.y,  false);
    }
}
void LineBuffer::write_distances_to(Polygon& subject, Polygon &other)
{
    for (Vertex v : subject.vertices())
    {
        int closest_edge;
        float closest_edge_alpha;
        float distance_from_other = distance(v.point, other, closest_edge, closest_edge_alpha);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << distance_from_other;
        g_renderer->get_font_renderer()->setColor(255, 255, 255);
        g_renderer->get_font_renderer()->addText(stream.str(), v.point.x, v.point.y,  false);
    }
}
void LineBuffer::append_lines_to_vector(Polygon& p)
{
    for (Edge e : p.edges()) {
        add_point(e.start);
        add_point(e.end);
    }
}
void LineBuffer::append_model_lines_to_vector(Polygon& p)
{
    for (Edge e : p.model_edges()) {
        add_point(e.start);
        add_point(e.end);
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
        add_point(vec_i);
        add_point(vec_j);
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
            add_point(point);
            add_point(next_point);
        } else { // Draw vector
            add_vector(point, next_point - point);
        }
    }

    /** Linear velocity **/
    add_vector(center, velocity);
}


void LineBuffer::add_point(glm::vec2 point)
{
    buffer.push_back(point.x);
    buffer.push_back(point.y);
    buffer.push_back(color.r);
    buffer.push_back(color.g);
    buffer.push_back(color.b);
}
void LineBuffer::add_point(float x, float y)
{
    buffer.push_back(x);
    buffer.push_back(y);
    buffer.push_back(color.r);
    buffer.push_back(color.g);
    buffer.push_back(color.b);
}


void LineBuffer::add_dot(glm::vec2 dot, float radius) 
{
    add_point(dot - glm::vec2(radius));
    add_point(dot + glm::vec2(radius));
    add_point(dot.x + radius, dot.y - radius);
    add_point(dot.x - radius, dot.y + radius);
}
void LineBuffer::add_dot(glm::vec2 dot)
{
    add_dot(dot, 2);
}
void LineBuffer::add_vector(glm::vec2 point, glm::vec2 vec)
{
    // const float radius = 0.07f;
    const float radius = glm::length(vec) / 15.f;
    const float arrow_angle = 2.7f;

    float vec_angle = atan2(vec.y, vec.x);
    glm::vec2 a1 = glm::vec2(cos(vec_angle - arrow_angle) * radius, sin(vec_angle - arrow_angle) * radius);
    glm::vec2 a2 = glm::vec2(cos(vec_angle + arrow_angle) * radius, sin(vec_angle + arrow_angle) * radius);

    add_point(point);
    add_point(point + vec);;

    add_point(point + vec);
    add_point(point + vec + a1);

    add_point(point + vec);
    add_point(point + vec + a2);
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
