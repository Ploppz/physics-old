#pragma once
#include "b.h"
class A
{
 public:
    A() {
        b = new B;
    };
    int c;
    B* b;
};
