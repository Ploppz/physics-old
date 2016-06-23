#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

class Debug;
typedef void (*manipulator)(Debug&);

// "doubt"

struct Beginning;

class Debug
{
 public:
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
    #define begin() _begin(__FUNCTION__, __FILE__, __LINE__);
    #define begin_c(condition) _begin(condition, __FUNCTION__, __FILE__, __LINE__);
    #define begin_d(description) _begin(description, __FUNCTION__, __FILE__, __LINE__);
    static void _begin(const std::string& function, const std::string& file, int line);
    static void _begin(bool condition, const std::string& function, const std::string& file, int line);
    static void _begin(const std::string& description, const std::string& function, const std::string& file, int line);
    static void end();
 private:

    static void prefix();
    static void draw_beginnings();
 public: /* Static members */
 private:
    static int level;
    static std::vector<Beginning> pending_beginnings;
    static bool next_is_newline;
    
    // implementation of conditional
    static bool do_not_write;
    static int do_not_write_counter;
};


template <>
inline Debug& Debug::operator <<(const char& rhs)
{
    // TODO Probably better to make this another function and let this overload be handled like the rest
    if (do_not_write) return *this;
    draw_beginnings();

    if (next_is_newline) {
        ++ level; prefix(); -- level;
        next_is_newline = false;
    }
    if (rhs == '\n') {
        next_is_newline = true;
    }
    std::cout << rhs;
    return *this;
}

// Manipulators
void newl(Debug& stream);
//


extern Debug dout;

