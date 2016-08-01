#pragma once
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include "PairManager.h"
#include "debug/debug.h"

template<int dimentions = 2>
struct AABB {
    AABB() : min {}, max {} {}
    float min[dimentions];
    float max[dimentions];
};
template <>
struct AABB<2> {
    AABB(float min_x, float max_x, float min_y, float max_y) {
        min[0] = min_x; min[1] = min_y;
        max[0] = max_x; max[1] = max_y;
    }
    AABB() : min{}, max{} {};
    float min[2];
    float max[2];
};

template <typename UserType, int dimentions = 2>
class SAP
{
 public: /* Structures */
    struct Box {
        int min_index[dimentions];
        int max_index[dimentions];
        UserType user_data;
    };
    struct EndPoint {
        EndPoint(int owner_index, float value, bool is_min) : owner_index(owner_index), value(value), is_min(is_min) {}
        int owner_index;
        float value;
        bool is_min;
    };
 public: /* Interface */
    SAP() : pairs {}, boxes {}, endpoints {} {};
    SAP(const std::vector<std::pair<AABB<dimentions>, UserType>>& boxes);
    UserType get_box_user_data(int box_index)
    {
        return boxes[box_index].user_data;
    }

    int add_box(AABB<dimentions> box, UserType user_data)
    {
        int new_box_index = boxes.size();
        for (int d = 0; d < dimentions; d ++) {
            float biggest_value = (endpoints[d].size() > 0) ? endpoints[d].back().value : 0;
            EndPoint min(new_box_index, biggest_value, true);
            EndPoint max(new_box_index, biggest_value, false);
            endpoints[d].push_back(min);
            endpoints[d].push_back(max);
        }

        /* Add new box and update it */
        // Thought: Maybe we can just insert the first axes where they should be.
        // The last axis, we place far away and then we update normally
        Box new_box;
        for (int d = 0; d < dimentions; d ++) {
             new_box.min_index[d] = endpoints[d].size() - 2; 
             new_box.max_index[d] = new_box.min_index[d] + 1; 
        }
        new_box.user_data = user_data;
        boxes.push_back(new_box);
        update_box(boxes.size() - 1, box);
        return new_box_index;
    }

    void update_box(int which, AABB<dimentions> new_box)
    {
        DebugBegin();
        for (int d = 0; d < dimentions; d ++)
        {
            Box& box = boxes[which];
            update_box_endpoint(d, box.min_index[d], new_box.min[d]);
            update_box_endpoint(d, box.max_index[d], new_box.max[d]);
        }
    }
    void remove_box(int which);

    void print()
    {
        for (int d = 0; d < dimentions; d ++)
        {
            std::cout << "Dimention " << d << ": " << std::endl;
            std::cout << " - Endpoints: ";
            for (auto it = endpoints[d].begin(); it != endpoints[d].end(); it ++)
            {
                std::cout << it->owner_index;
                if (it->is_min)
                    std::cout << "_ ";
                else
                    std::cout <<"^ ";
            }
            std::cout << std::endl;
        }
        for (int i = 0; i < boxes.size(); i ++ )
        {
            std::cout << i << ": x(" << endpoints[0][boxes[i].min_index[0]].value << ", " << endpoints[0][boxes[i].max_index[0]].value << "), y("
                << endpoints[1][boxes[i].min_index[1]].value << ", " << endpoints[1][boxes[i].max_index[1]].value << ")" << std::endl;
        }
    }
    void assert_validity()
    {
        for (int i = 0; i < boxes.size(); i ++ )
        {
            for (int d = 0; d < dimentions; d ++)
            {
                assert(endpoints[d][boxes[i].min_index[d]].is_min);
                assert(!endpoints[d][boxes[i].max_index[d]].is_min);
            }
        }
        for (int i = 0; i < boxes.size(); i ++ ) {
            for (int j = 0; j < boxes.size(); j ++ ) {
                if (i == j) continue;
                for (int d = 0; d < dimentions; d ++) {
                    assert(boxes[i].min_index[d] != boxes[j].min_index[d]);
                    assert(boxes[i].max_index[d] != boxes[j].max_index[d]);
                }  
            }
        }
    }

 private:
    int sign(float val) {
        return (0 < val) - (val < 0);
    }
    void update_box_endpoint(int d, int endpoint_index, float new_value)
    {
        DebugBeginC(false);
        EndPoint endpoint_copy = endpoints[d][endpoint_index]; // for insertion in the end
        float old_value = endpoint_copy.value;
        dout << "old value: " << old_value << newl;
        dout << "new value: " << new_value << newl;
        int direction = sign(new_value - old_value);
        dout << "direction: " << direction << newl;
        if (direction != 0) {
            int endpoint_new_index = endpoint_index + direction;
            while (endpoint_new_index >= 0 && endpoint_new_index < endpoints[d].size()) { // can be optimized - it's actually only for the initial round
                int new_direction = sign(new_value - endpoints[d][endpoint_new_index].value);
                dout << "value at next: " << endpoints[d][endpoint_new_index].value << newl;

                if (new_direction != direction)
                    break;

                // Signal that we jumped over this one
                if (endpoints[d][endpoint_new_index].is_min != endpoint_copy.is_min) {
                    // 'expands' is true if [a MIN endpoint moves left] or [a MAX endpoint moves right]
                    bool expands = (endpoint_copy.is_min) ^ (direction == 1);
                    if (expands)
                        overlap_starts_in_dimention(d, endpoint_copy.owner_index, endpoints[d][endpoint_new_index].owner_index);
                    else
                        overlap_ends_in_dimention(d, endpoint_copy.owner_index, endpoints[d][endpoint_new_index].owner_index);
                }
                    /* reevaluate_pair(endpoint_copy.owner_index, endpoints[d][endpoint_new_index].owner_index); */

                // Move this one behind us
                endpoints[d][endpoint_new_index - direction] = endpoints[d][endpoint_new_index];
                update_endpoint_index(d, endpoint_new_index - direction);
                dout << "endpoints[" << (endpoint_new_index - direction) << "] = " << endpoints[d][endpoint_new_index].value << newl;

                endpoint_new_index += direction;
                if (endpoint_new_index < 0 || endpoint_new_index >= endpoints[d].size())
                    break;
            }
            endpoint_new_index -= direction;
            endpoint_copy.value = new_value;
            endpoints[d][endpoint_new_index] = endpoint_copy;
            update_endpoint_index(d, endpoint_new_index);
            dout << "endpoints[" << endpoint_new_index << "] = " << endpoint_copy.value << newl;
            dout << "final index: " << endpoint_new_index << newl;

            
            /* Insert */
            // move_to(endpoint_index, endpoint_new_index, endpoints[d]);
        }
    }
    /* Checks all _other_ dimentions for overlap*/
    void overlap_starts_in_dimention(int dimention, int box1_index, int box2_index)
    {
        DebugBeginC(false);
        bool overlap = true;
        for (int d = 0; d < dimentions; d ++)
        {
            if (d != dimention) {
                overlap = overlap && overlap_intervals( boxes[box1_index].min_index[d],  boxes[box1_index].max_index[d],
                                                    boxes[box2_index].min_index[d],  boxes[box2_index].max_index[d]);
                if (overlap) dout << " YES overlap: ";
                else        dout << " NO overlap: ";
                dout << "(" << boxes[box1_index].min_index[d] << ", " << boxes[box1_index].max_index[d] << ") vs ("
                    << boxes[box2_index].min_index[d] << ", " << boxes[box2_index].max_index[d] << ")" << newl;
            }
        }
        if (overlap)
            pairs.register_pair(box1_index, box2_index);
        dout << "State of x: " << newl;
        for (EndPoint p : endpoints[0]) {
            dout << p.value << ", ";
        }
    }

    // 'the dimention isn't really necessary here - we could just always ask PairManager to remove it. But this is maybe more performant?
    void overlap_ends_in_dimention(int dimention, int box1_index, int box2_index)
    {
        pairs.unregister_pair(box1_index, box2_index);
    }

    bool overlap_intervals(int box1_min, int box1_max, int box2_min, int box2_max)
    {
        return (box1_min > box2_min && box1_min < box2_max) || (box2_min > box1_min && box2_min < box1_max);
        // TODO <= and >= - however, also fix so that boxes aren't tested with themselves...
        // or don't really need = since indices aren't supposed to be the same...
    }

    void update_endpoint_index(int d, int endpoint_index) // Writes to the owner
    {
        if (endpoints[d][endpoint_index].is_min)
            boxes[  endpoints[d][endpoint_index].owner_index  ].min_index[d] = endpoint_index;
        else
            boxes[  endpoints[d][endpoint_index].owner_index  ].max_index[d] = endpoint_index;
    }
    

    void move_to(int src, int dest, std::vector<EndPoint>& container) /* UNUSED */
    {
        EndPoint tmp = container[src];

        int direction = sign(dest - src);
        for (int i = src + direction; i != dest; i += direction)
        {
            container[i] = container[i - direction];
        }
        container[dest] = tmp;
    }


 public: /* Data */
    PairManager pairs;
 private: /* Data */
    std::vector<Box> boxes;
    std::vector<EndPoint> endpoints[dimentions];
};
