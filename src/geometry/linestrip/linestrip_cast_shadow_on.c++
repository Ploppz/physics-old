#include "LineStrip.h"
#include "../Intersection.h"
#include "../geometry.h"
#include <limits>
#include <glm/glm.hpp>
using namespace glm;

LineStrip Intersection::cast_shadow_on(Polygon* polygon, glm::vec2 direction)
{
    // Align transform makes `direction` toward +x
    mat2 align_transform = rotate_coor_system(direction);
            mat2 align_inverse = glm::inverse(align_transform);
    // For either of y_min and y_max vertices: If [vertex is intersect]:
    //   use that vertex
    // else:
    //   Raycast from the vertex, find the two intersects with minimum x (closest to Intersection)
    //    - if no intersect: use the polygon's y_min and/or y_max vertices
    //    (- if found none of the two intersect: test if the polygon is at all between the two rays for debug)
    // 
    // LineStrip extraction:
    //   at the start, the polygon needs to go toward the shadow.. that is, geometrically toward the other vertex
    int this_min_y_vertex,
        this_max_y_vertex;
    float this_min_y =      std::numeric_limits<float>::max(),
          this_max_y =    - std::numeric_limits<float>::max();
    bool using_this_min_y_vertex = false,
         using_this_max_y_vertex = false;

    for (int i = 0; i < vertices.size(); i ++)
    {
        vec2 transformed_vertex = align_transform * vertices[i].point;
        if (transformed_vertex.y < this_min_y) {
            this_min_y = transformed_vertex.y;
            this_min_y_vertex = i;
            if (vertices[i].intersect)
                using_this_min_y_vertex = true;
        }
        if (transformed_vertex.y > this_max_y) {
            this_max_y = transformed_vertex.y;
            this_max_y_vertex = i;
            if (vertices[i].intersect)
                using_this_max_y_vertex = true;
        }
    }

    int polygon_min_y_min_x_vertex,
        polygon_max_y_min_x_vertex;
    float polygon_min_y_min_x_alpha,
          polygon_max_y_min_x_alpha;
    float polygon_min_y_min_x = std::numeric_limits<float>::max(),
          polygon_max_y_min_x = std::numeric_limits<float>::max();
    Polygon::Edge polygon_edge_it = polygon->first_edge(),
                  polygon_end_it  = polygon->last_edge();

    while (true) {
        vec2 transformed_start = align_transform * polygon_edge_it.start_tr();
        vec2 transformed_end   = align_transform * polygon_edge_it.end_tr();
        /* min_y_vertex ray */
        if ( ! using_this_min_y_vertex) {
            float alpha;
            float intersect_x =
                intersect_horizontal(transformed_start, transformed_end - transformed_start, this_min_y, alpha);
            if (alpha >= 0 && alpha < 1) {
                if (intersect_x < polygon_min_y_min_x) {
                    polygon_min_y_min_x_vertex = polygon_edge_it.get_index();
                    polygon_min_y_min_x_alpha = alpha;
                    polygon_min_y_min_x = intersect_x;
                }
            }
        }
        /* max_y_vertex ray */
        if ( ! using_this_max_y_vertex) {
            float alpha;
            float intersect_x =
                intersect_horizontal(transformed_start, transformed_end - transformed_start, this_max_y, alpha);
            if (alpha >= 0 && alpha < 1) {
                if (intersect_x < polygon_max_y_min_x) {
                    polygon_max_y_min_x_vertex = polygon_edge_it.get_index();
                    polygon_max_y_min_x_alpha = alpha;
                    polygon_max_y_min_x = intersect_x;
                }
            }
        }
        /** Iteration **/
        if (polygon_edge_it == polygon_end_it)
            break;
        ++ polygon_edge_it;
    }
    if (using_this_min_y_vertex) {
        const HybridVertex& vertex = vertices[this_min_y_vertex];
        if (vertex.edge1_owner == polygon) {
            polygon_min_y_min_x_vertex = vertex.edge1_index;
            polygon_min_y_min_x_alpha = vertex.alpha1;
        } else {
            polygon_min_y_min_x_vertex = vertex.edge2_index;
            polygon_min_y_min_x_alpha = vertex.alpha2;
        }
    }
    if (using_this_max_y_vertex) {
        const HybridVertex& vertex = vertices[this_max_y_vertex];
        if (vertex.edge1_owner == polygon) {
            polygon_max_y_min_x_vertex = vertex.edge1_index;
            polygon_max_y_min_x_alpha = vertex.alpha1;
        } else {
            polygon_max_y_min_x_vertex = vertex.edge2_index;
            polygon_max_y_min_x_alpha = vertex.alpha2;
        }
    }

    bool min_y_suitable_for_start = false;
    // min_y vertex is suitable for start of linestrip if the direction if the transformed
    // polygon at that point goes toward 
    // or... If the polygon at min_y_vertex goes toward positive y!
    vec2 polygon_min_y_min_x_tr =      align_transform * polygon->transformed(polygon_min_y_min_x_vertex);
    vec2 polygon_min_y_min_x_next_tr = align_transform * polygon->transformed(polygon_min_y_min_x_vertex + 1);

    if ( (polygon_min_y_min_x_next_tr - polygon_min_y_min_x_tr).x > 0) /* edge goes upward */
    {
        min_y_suitable_for_start = true;
    }
    
    LineStrip result(polygon);
    result.set_start<false>( EdgePoint(polygon_min_y_min_x_vertex, polygon_min_y_min_x_alpha, polygon) );
    result.set_start<false>( EdgePoint(polygon_max_y_min_x_vertex, polygon_max_y_min_x_alpha, polygon) );
    if ( ! min_y_suitable_for_start) {
        result.swap_start_and_end();
    }
    return result;
}
