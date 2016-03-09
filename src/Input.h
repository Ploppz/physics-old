#include <GLFW/glfw3.h>


class Input
{
public:
    static void UpdateMouse(GLFWwindow* window);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

    /* Variables */
    static bool keys[512];
    static bool mouse_left;
    static int mouse_x;
    static int mouse_y;
    static double scroll;

    static int mouse_drag_x;
    static int mouse_drag_y;
    
private:

    static int prev_mouse_x;
    static int prev_mouse_y;
    // No dynamic class
    Input();
};
