#include "Statistics.h"

#include <ostream>

#define USING_EXP_AVG true

Statistics::Statistics(float initial_value)
{
    add_value(initial_value);
}
void Statistics::add_value(int value)
{
    add_value((float)value);
}
void Statistics::add_value(float value)
{
    if (num_samples == 0) {
        avg = value;
        ++ num_samples;
    } else {
        if (USING_EXP_AVG) {
            add_exp_avg(value, avg, 0.2f);
        } else {
            add_avg(value, avg, num_samples);
        }
    }

    if (value > max)
        max = value;
    if (value < min)
        min = value;
    if (history.size() > history_length)
    {
        history.pop_front();
    }
    history.push_back(value);
}


////
std::ostream& operator<< (std::ostream& out, Statistics& s)
{
    out << "{avg=" << s.get_avg() << ", max=" << s.get_max() << ", min=" << s.get_min() << "}";
    return out;
}
