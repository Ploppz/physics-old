#pragma once
#include <list>
#include "LineStrip.h"
#include "LineStripSeries.h"

template <bool reverse>
struct linestrip_iterator;

template <>
struct linestrip_iterator<false>
{
    std::list<LineStrip>::iterator iterator;
//
    std::list<LineStrip>::iterator base() {
        return iterator;   
    }
    void insert_after_and_move_to_new(LineStrip linestrip, std::list<LineStrip>& container) {
        ++ iterator;
        container.insert( iterator, linestrip);
        -- iterator;
    }
    void set_to_start_of(std::list<LineStrip>& container) {
        iterator = container.begin();   
    }
    void set_to_end_of(std::list<LineStrip>& container) {
        iterator = container.end();   
    }
    bool at_start_of(std::list<LineStrip>& container) {
        return iterator == container.begin();
    }

    bool at_end_of(std::list<LineStrip>& container) {
        return iterator == (--container.end());
    }

};
template <>
struct linestrip_iterator<true>
{

    std::list<LineStrip>::reverse_iterator iterator;
//
    std::list<LineStrip>::iterator base() {
        return  -- iterator.base();   
    }

    void insert_after_and_move_to_new(LineStrip linestrip, std::list<LineStrip>& container) {
        container.insert( -- iterator.base(), linestrip);
        ++ iterator;
    }

    void set_to_start_of(std::list<LineStrip>& container) {
        iterator = container.rbegin();   
    }
    void set_to_end_of(std::list<LineStrip>& container) {
        iterator = container.rend();   
    }
    bool at_start_of(std::list<LineStrip>& container) {
        return iterator == container.rbegin();
    }

    bool at_end_of(std::list<LineStrip>& container) {
        return iterator == (--container.rend());
    }
};
