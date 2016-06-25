#pragma once
class Body;

class ResolutionAlg
{
public:
    virtual void init() = 0;
    virtual void iteration_start() = 0;
    virtual void treat(Body, Body, float delta_time) = 0;
    virtual bool done() = 0;
};
