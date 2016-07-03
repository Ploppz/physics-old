#pragma once
#include "Statistics.h"

#include <map>

// Initially, I had only 'statistics_objects'. Don't know what the best
// way to also have 'counters' is
class StatisticsCollection
{
 public:
    StatisticsCollection() {}
    void add_value(const std::string& name, float value);
    void count(const std::string& name);
    /* Iterating */
    std::map<std::string, Statistics>::iterator begin();
    std::map<std::string, Statistics>::iterator end();
    std::map<std::string, int>::iterator begin_counters();
    std::map<std::string, int>::iterator end_counters();
 private:
    std::map<std::string, Statistics> statistics_objects;
    std::map<std::string, int> counters;
};
