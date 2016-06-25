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

void a()
{
    DebugBegin();
    dout << "a called" << newl;
}
void b()
{
    DebugBegin();
    dout << "b called" << newl;
    a();
}
void c()
{
    DebugBegin();
    dout << "c called" << newl;
    b();
    b();
}
void lala()
{
    DebugBegin();
    dout << "lala called" << newl;
    c();
    dout << "lala: test" << newl;
    b();
    a();
}

int main()
{
    Debug::set_color1(brown);
    Debug::set_color2(gray);

    DebugBegin();
    dout << "one test bla" << newl;
    lala();
    A a;
    a.x = 3;
    dout << a << newl;
}
