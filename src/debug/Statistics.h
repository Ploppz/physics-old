#pragma once
#include <deque>
#include <limits>
#include <ostream>

class Statistics
{
 public:
    Statistics() {};
    Statistics(float initial_value);
    void add_value(float value);
    void add_value(int value);

    /* Getters */
    float get_avg() { return avg; }
    float get_max() { return max; }
    float get_min() { return min; }
    const std::deque<float>& get_history() { return history; }
 private:
    void add_avg(float new_value, float& existing_value, int& num_samples)
    {
        existing_value = (new_value + existing_value * num_samples) / (num_samples + 1);
        ++ num_samples;
    }

    void add_exp_avg(float new_value, float& existing_value, float weight_of_new)
    {
        existing_value = weight_of_new * new_value + (1 - weight_of_new) * existing_value;
        num_samples = 1;
    }
 private: /* Members */
    float avg = 0;
    int num_samples = 0;
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::min();
    std::deque<float> history;
    unsigned int history_length = 20;
};

std::ostream& operator<< (std::ostream& out, Statistics& s);
