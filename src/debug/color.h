#pragma once
#include <iostream>

typedef std::ostream& (*std_manipulator)(std::ostream&);

//TODO next cannot declare it using typedef
std::ostream& nocolor(std::ostream&);

std::ostream& red(std::ostream&);
std::ostream& green(std::ostream&);
std::ostream& brown(std::ostream&);
std::ostream& blue(std::ostream&);
std::ostream& magenta(std::ostream&);
std::ostream& cyan(std::ostream&);
std::ostream& gray(std::ostream&);

std::ostream& Red(std::ostream&);
std::ostream& Green(std::ostream&);
std::ostream& Brown(std::ostream&);
std::ostream& Blue(std::ostream&);
std::ostream& Magenta(std::ostream&);
std::ostream& Cyan(std::ostream&);
std::ostream& Gray(std::ostream&);
