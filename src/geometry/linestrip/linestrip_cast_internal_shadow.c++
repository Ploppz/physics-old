#include "LineStrip.h"
#include "../Intersection.h"
#include <glm/glm.hpp>

using namespace glm;
LineStrip Intersection::cast_internal_shadow(glm::vec2 direction, Polygon* subject, Renderer &renderer)
{
    direction = normalize(direction);
    float align_transform_d[4] = {direction.x, -direction.y, direction.y, direction.x}; //note: column major
    mat2 align_transform = make_mat2x2(align_transform_d);
    mat2 align_inverse = glm::inverse(align_transform);

    /* find the limits of the intersection on the transformed y-axis */
    float y_min = std::numeric_limits<float>::max();
    float y_min_x;
    float y_max = -std::numeric_limits<float>::max();
    Vertex y_min_vertex, y_max_vertex;

    for (int i = 0; i < vertices.size(); i ++) 
    {
        vec2 vertex = align_transform * vertices[i].point;
        if (vertex.y < y_min) {
            y_min = vertex.y;
            y_min_x = vertex.x;
            y_min_vertex = Vertex(i, this);
        }
        if (vertex.y > y_max) {
            y_max = vertex.y;
            y_max_vertex = Vertex(i, this);
        }
    }
    // renderer.add_dot(align_inverse * vec2(y_min_x, y_max));

    float alpha;

    EdgePoint start_point(subject);
    float start_x = - 1; // for finding the smallest x value > 0
    EdgePoint end_point(subject);
    bool end_found = false;
    float end_best_x = std::numeric_limits<float>::max();

    // the following is needed because we continously look for improvement, but we are not sure
    // whether we will find an end once we have found a start
    EdgePoint start_point_wip(subject);
    bool start_found_wip = false;


    // (i, j) are the edge indeces, (a, b) are the edge coordinates
    int i = subject->vertices.size() - 1; 
    vec2 a = align_transform * subject->transform(subject->vertices.back());

    for (int j_unbounded = 0; j_unbounded < subject->vertices.size() * 2; j_unbounded ++) /*[1]*/
    {
        int j = j_unbounded % subject->vertices.size();
        vec2 b = align_transform * subject->transformed(j);
        
        if (a.y >= y_max && b.y <= y_max)
        {
            // Enter via y_max - start recording the path
            intersect_horizontal(a, b - a, y_max, alpha);
            start_found_wip = true;
            start_point_wip.index = i;
            start_point_wip.alpha = alpha;
        }
        if (a.y >= y_min && b.y <= y_min)
        {
            // Leave via y_min
            if (start_found_wip)
            {
                intersect_horizontal(a, b - a, y_min, alpha);
                float x = (a * (1 - alpha) + b * alpha).x;
                if (x >= y_min_x - 1 && x < end_best_x) /*[2]*/
                {
                    start_point.index = start_point_wip.index;
                    start_point.alpha = start_point_wip.alpha;
                    end_point.index = i;
                    end_point.alpha = alpha;

                    end_best_x = x;
                    end_found = true;
                    start_found_wip = false;
                }
            }
        }
        i = j;
        a = b;
    }
    assert(end_found);
    /* Create line strip */
	bool CCW = (subject->signed_area() > 0); // Is this intersection counter clockwise?
    LineStrip result(subject);
    // result.CCW = CCW;
    if (CCW) {
        result.set_start<true>(start_point);
        result.set_end<true>(end_point);
    } else {
        result.set_start<false>(start_point);
        result.set_end<false>(end_point);
    }
    return result;
}
