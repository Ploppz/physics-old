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



Debug dout;
int                     Debug::level                = -1;
bool                    Debug::next_is_newline      = true;
std::vector<Beginning>  Debug::pending_beginnings;
bool Debug::do_not_write = false;
int  Debug::do_not_write_counter = 0;



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
    for (Beginning& b : pending_beginnings)
    {
        ++ level;
        prefix();
        std::cout << "┌─── " << b.function  << "() ──── " << b.file << ":" << b.line << std::endl;
        next_is_newline = true;
    }
    pending_beginnings.clear();
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
    if (!next_is_newline)
        std::cout << std::endl;
    next_is_newline = true;
    prefix();
    std::cout << "└──────────────────────────────────────────────────" << std::endl; 

    -- level;
    assert(level >= -1);
}
void Debug::prefix()
{
    for (int i = 0; i < level; i ++)
    {
        std::cout << "│  ";
    }
}

Debug& Debug::operator << (const manipulator m)
{
    (*m)(*this);
    return *this;
}

void newl(Debug& stream)
{
    stream << "\n";
}
