#include <iostream>
#include "../src/debug.h"
#include "../src/render/Renderer.h"

FontRenderer *g_font_renderer;
Renderer *g_renderer;

struct A
{
    int x;
};
std::ostream& operator << (std::ostream& lhs, A rhs)
{
    lhs << "A.x = " << rhs.x;
    return lhs;
};

void lala()
{
    Debug::begin();
    int a = 3;
    dout << "lala test" << a << " 2 .. " << 1 << newl;
    Debug::begin();
    Debug::begin();
    dout << "hei";
    Debug::end();
    Debug::begin();
    dout << "hoi";
    Debug::end();
    Debug::begin();
    dout << "hi";
    Debug::end();
    Debug::end();

    Debug::end();
}

int main()
{
    Debug::begin();
    dout << "one test bla" << newl;
    lala();
    A a;
    a.x = 3;
    dout << a << newl;
    Debug::end();
}
