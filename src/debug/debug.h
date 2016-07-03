#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "color.h"

class Debug;
typedef void (*manipulator)(Debug&);

/** Creates a Debug object ... RAII style **/
#define DebugBegin() Debug dout(__FUNCTION__, __FILE__, __LINE__);
#define DebugBeginC(condition) Debug dout(condition, __FUNCTION__, __FILE__, __LINE__);

// "doubt"

struct Beginning;

class Debug
{
 public:
    Debug(const std::string& function, const std::string& file, int line);
    Debug(bool condition, const std::string& function, const std::string& file, int line);
    ~Debug();

    void fatal(const std::string& message);

    Debug& operator << (const manipulator);

    template <typename T>
    Debug& operator << (const T& rhs)
    {
        std::ostringstream ss;
        ss << rhs;
        std::string s = ss.str();
        for (char c : s)
            *this << c;
        return *this;
    }


 public: /* Static methods */
    static void _begin(const std::string& function, const std::string& file, int line);
    static void _begin(bool condition, const std::string& function, const std::string& file, int line);
    static void _begin(const std::string& description, const std::string& function, const std::string& file, int line);
    static void end();

    static void set_color1(std_manipulator color);
    static void set_color2(std_manipulator color);
 private:

    static void prefix();
    static void prefix(int level);
    static void draw_beginnings();
    static void draw_end();
 public: /* Static members */
 private:
    static int level;
    static std::vector<Beginning> pending_beginnings;
    static bool next_is_newline;
    
    // sibling detection
    static int pending_ends;
    // implementation of conditional
    static bool do_not_write;
    static int do_not_write_counter;
    // config..
    static std_manipulator color1;
    static std_manipulator color2;
};


template <>
inline Debug& Debug::operator <<(const char& rhs)
{
    // TODO Probably better to make this another function and let this overload be handled like the rest
    if (do_not_write) return *this;
    // only one of these two functions will do something.
    draw_beginnings();
    draw_end();

    /* std::cout << "$";  */
    if (next_is_newline) {
        ++ level; prefix(); -- level;
        next_is_newline = false;
    }
    if (rhs == '\n') {
        next_is_newline = true;
    }
    /* std::cout << "(" << rhs << ")"; */
    std::cout << rhs;
    return *this;
}

// Manipulators
void newl(Debug& stream);
//


extern Debug dout;

