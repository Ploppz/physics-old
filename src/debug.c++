#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <limits>

// all overloaded <<, write to string stream, then iteratively call the overload for char
// TODO: don't write empty boxes..
// ... will need some sort of FIFO
// might not need wrote_since_beginning

#include "debug.h"
#include "color.h"

/*** Dynamic ***/
Debug::Debug(const std::string& function, const std::string& file, int line) {
    Debug::_begin(function, file, line);
}
Debug::Debug(bool cond, const std::string& function, const std::string& file, int line) {
    Debug::_begin(cond, function, file, line);
}
Debug::~Debug() {
    Debug::end();
}


/*** Static ***/

struct Beginning
{
    Beginning(const std::string& function, const std::string& file, int line)
        : function(function), file(file), line(line) {}
    Beginning(const std::string& function, const std::string& file, int line, bool is_sibling)
        : function(function), file(file), line(line), is_sibling(is_sibling) {}
    std::string function;
    std::string file;
    int line;
    bool is_sibling = false;
};



int                     Debug::level                = -1;
bool                    Debug::next_is_newline      = true;
//
std::vector<Beginning>  Debug::pending_beginnings;
int Debug::pending_ends = 0;
//
bool Debug::do_not_write = false;
int  Debug::do_not_write_counter = 0;
//
std_manipulator Debug::color1 = brown;
std_manipulator Debug::color2 = gray;



void Debug::_begin(const std::string& function, const std::string& file, int line)
{
    if (do_not_write) {
        ++ do_not_write_counter;
    } else {
        pending_beginnings.push_back(Beginning(function, file, line));
    }
}
void Debug::_begin(bool cond, const std::string& function, const std::string& file, int line)
{
    if (!cond) {
        if (do_not_write) {
            ++ do_not_write_counter;
        } else {
            do_not_write_counter = 1;
            do_not_write = true;
        }
    } else {
        pending_beginnings.push_back(Beginning(function, file, line));
    }
}
void Debug::draw_beginnings()
{
   if (pending_beginnings.size() == 0) return;
    if (!next_is_newline)
        std::cout << std::endl;
    next_is_newline = true;
   if (pending_ends == 1) {
       // TODO wrong order?? need to get from front?
        Beginning b = pending_beginnings.back();
        pending_beginnings.pop_back();
        next_is_newline = true;
        prefix();
        std::cout << color2 << "├─── "
                  << color1 << b.function
                  << color2 << "() ──── "
                  << color1 << b.file
                  << color2 << ":"
                  << color1 << b.line
                  << nocolor << std::endl;
        -- pending_ends;
    }

    draw_end();
    for (Beginning& b : pending_beginnings)
    {
        ++ level;
        prefix();
        std::cout << color2 << "┌─── "
                  << color1 << b.function
                  << color2 << "() ──── "
                  << color1 << b.file
                  << color2 << ":"
                  << color1 << b.line
                  << nocolor << std::endl;
        next_is_newline = true;
    }
    pending_beginnings.clear();
}
void Debug::draw_end()
{
    if (pending_ends == 0) return;
    if (!next_is_newline)
        std::cout << std::endl;
    next_is_newline = true;

    int start_level = level - pending_ends + 1;

    prefix(start_level);
    std::cout << color2; // TODO needed?

    { /* First end */
        -- level;
        -- pending_ends;
        std::cout << "└──";
    }
    while (pending_ends > 0)
    {
        -- level;
        -- pending_ends;
        std::cout << "┴──";
    }
    std::cout << "──────────────────────────────" << std::endl;
    pending_ends = 0;

    assert(level >= -1);
}

void Debug::end()
{
    if (do_not_write) {
        -- do_not_write_counter;
        if (do_not_write_counter == 0) {
            do_not_write = false;
        }
        return;
    }
    if (pending_beginnings.size() > 0) {
        pending_beginnings.pop_back();
        return;
    }
    /* Draw end.. or that is, make it happen some time */
    ++ pending_ends;
}
void Debug::prefix()
{
    std::cout << color2;
    for (int i = 0; i < level; i ++)
    {
        std::cout << "│  ";
    }
    std::cout << nocolor;
}
void Debug::prefix(int level)
{
    std::cout << color2;
    for (int i = 0; i < level; i ++)
    {
        std::cout << "│  ";
    }
    std::cout << nocolor;
}

void Debug::set_color1(std_manipulator color) { color1 = color; }
void Debug::set_color2(std_manipulator color) { color2 = color; }
Debug& Debug::operator << (const manipulator m)
{
    (*m)(*this);
    return *this;
}

void newl(Debug& stream)
{
    stream << "\n";
}
