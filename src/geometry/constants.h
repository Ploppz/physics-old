#pragma once
enum Axis { X=0, Y=1 };


typedef bool Side;
const bool OUT = false; // = BACK
const bool IN = true; // = FORTH

enum Direction {
    BACK = false,  FORTH = true, // previously LEFT and RIGHT
    LEFT = BACK, RIGHT = FORTH
};

