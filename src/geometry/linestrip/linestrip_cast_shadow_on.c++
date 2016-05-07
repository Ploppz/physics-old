#include "LineStrip.h"
#include "../Intersection.h"
#include "../geometry.h"
#include <limits>
#include <glm/glm.hpp>
using namespace glm;

LineStrip Intersection::cast_shadow_on(Polygon* polygon, glm::vec2 direction)
{
    const bool DEBUG = false;
    // Align transform makes `direction` toward +x
    mat2 align_transform = rotate_coor_system(direction);
    mat2 align_inverse = glm::inverse(align_transform);
    // For either of y_min and y_max vertices: If [vertex includes *polygon]:
    //   use that vertex
    // else:
    //   Raycast from the vertex, find the two intersects with minimum x > [x of ray origin]
    //    - if no intersect: use the polygon's y_min and/or y_max vertices
    //    (- if found none of the two intersect: test if the polygon is at all between the two rays for debug)
    // 
    // LineStrip extraction:
    //   at the start, the polygon needs to go toward the shadow.. that is, geometrically toward the other vertex
    

    /*** Find min_y and max_y vertices of Intersection ***/

    int this_min_y_index,
        this_max_y_index;
    /* The minimum and maximum-y vertices of the Intersection, transformed by align transform */
    glm::vec2 this_min_y_vec(0, std::numeric_limits<float>::max()),
              this_max_y_vec(0, - std::numeric_limits<float>::max());
    bool using_this_min_y_vertex = false,
         using_this_max_y_vertex = false;

    for (int i = 0; i < vertices.size(); i ++)
    {
        vec2 transformed_vertex = align_transform * vertices[i].point;
        if (transformed_vertex.y < this_min_y_vec.y) {
            /* UPDATE min_y vertex */
            this_min_y_vec = transformed_vertex;
            this_min_y_index = i;
            /* this_min_y_x = transformed_vertex.y */
            if (vertices[i].intersect || vertices[i].owner == polygon) { // if *polygon is included in vertex
                using_this_min_y_vertex = true;
            } else {
                using_this_min_y_vertex = false;
            }

        }
        if (transformed_vertex.y > this_max_y_vec.y) {
            /* UPDATE max_y vertex */
            this_max_y_vec = transformed_vertex;
            this_max_y_index = i;
            if (vertices[i].intersect || vertices[i].owner == polygon) { // if *polygon is included in vertex
                using_this_max_y_vertex = true;
            } else {
                using_this_max_y_vertex = false;
            }
        }
    }
    if (DEBUG) std::cout << "Using this: " << std::boolalpha << using_this_min_y_vertex << ", " << using_this_max_y_vertex << std::endl;

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
                intersect_horizontal(transformed_start, transformed_end - transformed_start, this_min_y_vec.y, alpha);
            if (intersect_x >= this_min_y_vec.x) { /* right of the start of ray */
                if (alpha >= 0 && alpha < 1) {
                    if (DEBUG) std::cout << "Intersect (min) ";
                    if (DEBUG) std::cout << " (" << intersect_x << ") ";
                    if (intersect_x < polygon_min_y_min_x) {
                        if (DEBUG) std::cout << "got accepted" << std::endl;
                        polygon_min_y_min_x_vertex = polygon_edge_it.get_index();
                        polygon_min_y_min_x_alpha = alpha;
                        polygon_min_y_min_x = intersect_x;
                    }
                    if (DEBUG) std::cout << std::endl;
                }
            }
        }
        /* max_y_vertex ray */
        if ( ! using_this_max_y_vertex) {
            float alpha;
            float intersect_x =
                intersect_horizontal(transformed_start, transformed_end - transformed_start, this_max_y_vec.y, alpha);
            if (intersect_x >= this_max_y_vec.x) { /* right of the start of ray */
                if (alpha >= 0 && alpha < 1) {
                    if (DEBUG) std::cout << "Intersect (max) ";
                    if (DEBUG) std::cout << " (" << intersect_x << ") ";
                    if (intersect_x < polygon_max_y_min_x) {
                        if (DEBUG) std::cout << "got accepted" << std::endl;
                        polygon_max_y_min_x_vertex = polygon_edge_it.get_index();
                        polygon_max_y_min_x_alpha = alpha;
                        polygon_max_y_min_x = intersect_x;
                    }
                    if (DEBUG) std::cout << std::endl;
                }
            }
        }
        /** Iteration **/
        if (polygon_edge_it == polygon_end_it)
            break;
        ++ polygon_edge_it;
    }
    if (using_this_min_y_vertex) {
        const HybridVertex& vertex = vertices[this_min_y_index];
        if (vertex.edge1_owner == polygon) {
            polygon_min_y_min_x_vertex = vertex.edge1_index;
            polygon_min_y_min_x_alpha = vertex.alpha1;
        } else {
            polygon_min_y_min_x_vertex = vertex.edge2_index;
            polygon_min_y_min_x_alpha = vertex.alpha2;
        }
    }
    if (using_this_max_y_vertex) {
        const HybridVertex& vertex = vertices[this_max_y_index];
        if (vertex.edge1_owner == polygon) {
            polygon_max_y_min_x_vertex = vertex.edge1_index;
            polygon_max_y_min_x_alpha = vertex.alpha1;
        } else {
            polygon_max_y_min_x_vertex = vertex.edge2_index;
            polygon_max_y_min_x_alpha = vertex.alpha2;
        }
    }

    // FLAW in following: min_y may be a local minima.
    bool min_y_suitable_for_start = false;
    // min_y vertex is suitable for start of linestrip if the direction if the transformed
    // polygon at that point goes toward 
    // or... If the polygon at min_y_vertex goes toward positive y!
    if ( ! (using_this_min_y_vertex && !vertices[this_min_y_index].intersect) ) {
        vec2 polygon_min_y_min_x_tr =      align_transform * polygon->transformed(polygon_min_y_min_x_vertex);
        vec2 polygon_min_y_min_x_next_tr = align_transform * polygon->transformed(polygon_min_y_min_x_vertex + 1);

        if ( (polygon_min_y_min_x_next_tr - polygon_min_y_min_x_tr).y > 0) /* edge goes upward */
        {
            min_y_suitable_for_start = true;
        }
    } else {
        // local minima: it's suitable if the next edge on the polygon points
        // more to the right (+x) than previous
        Polygon::Edge next_edge(polygon_min_y_min_x_vertex, polygon);
        Polygon::Edge prev_edge(polygon_min_y_min_x_vertex - 1, polygon);
        glm::vec2 next_edge_vec = (next_edge.end_tr() - next_edge.start_tr());
        glm::vec2 prev_edge_vec = (prev_edge.end_tr() - prev_edge.start_tr());
        min_y_suitable_for_start = !leftof(prev_edge_vec, next_edge_vec);
        if (DEBUG) std::cout << "using this min_y vertex" << std::endl;
    }
    
    LineStrip result(polygon);
    result.set_start<false>( EdgePoint(polygon_min_y_min_x_vertex, polygon_min_y_min_x_alpha, polygon) );
    result.set_end<false>( EdgePoint(polygon_max_y_min_x_vertex, polygon_max_y_min_x_alpha, polygon) );
    if ( ! min_y_suitable_for_start) {
        if (DEBUG) std::cout << "SWAP" << std::endl;
        result.swap_start_and_end();
    }
    if (DEBUG) std::cout << result << std::endl;
    return result;
}
