#pragma once
#include "ResolutionAlg.h"
#include "Contact.h"

class DepthResolutionAlg
{
 public:
    void init();
    void iteration_start();
    void treat(Body, Body, float delta_time);
    bool done();
 private:
    void resolve(Body, Body, DepthContact);
 private:
    float biggest_depth;
    int iterations;
 private:
    const float CORRECTION_STRENGTH = 0.7f;
    const float DEPTH_OK_THRESHOLD = 10;
    const int MAX_ITERATIONS = 12;
};
