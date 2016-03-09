#include "Input.h"
#include <algorithm>
#include <iostream>
#include <GLFW/glfw3.h>

bool Input::keys[512] {};
bool Input::mouse_left {};
int Input::mouse_x {};
int Input::mouse_y {};
double Input::scroll {};
int Input::mouse_drag_x {};
int Input::mouse_drag_y {};
int Input::prev_mouse_x {};
int Input::prev_mouse_y {};


void Input::UpdateMouse(GLFWwindow* window)
{
    /// mouse position
    double xpos, ypos;
    glfwGetCursorPos( window, &xpos, &ypos);
    mouse_x = xpos;
    mouse_y = ypos;

    /// mouse drag
    if (mouse_left) {
        mouse_drag_x = mouse_x - prev_mouse_x;
        mouse_drag_y = -(mouse_y - prev_mouse_y);
    } else {
        mouse_drag_x = mouse_drag_y = 0;
    }
    prev_mouse_x = mouse_x;
    prev_mouse_y = mouse_y;
    
}


void Input::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    scroll += yoffset;
}
void Input::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        keys[key] = true;
    }
    if (action == GLFW_RELEASE) {
        keys[key] = false;
    }
}
void Input::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouse_left = true;
        } else if (action == GLFW_RELEASE) {
            mouse_left = false;
        }
    }
}

