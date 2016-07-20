#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "debug/debug.h"
// Interface / abstract class for classes that can render things


class Renderer
{
 public:
    virtual void render(float center_x, float center_y, int width, int height, float zoom) = 0;

    void set_active_buffer( const std::string& name )
    {
        active_buffer = &buffers[name];
        active_buffer_name = name;
    }

    std::string get_active_buffer_name() { return active_buffer_name; }

    void use_default_buffer()
    {
        set_active_buffer("default");
    }
    void clear_buffer()
    {
        active_buffer->clear();
    }
 protected:
    std::unordered_map<std::string, std::vector<float>> buffers;
    std::vector<float>* active_buffer;
    std::string active_buffer_name; // primarily useful for debugging
};
