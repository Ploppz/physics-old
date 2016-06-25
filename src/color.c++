#include "color.h"
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

string color(int num)
{
    stringstream ss;
    ss << "\e[" << num << "m";
    return ss.str();
}
std::ostream& nocolor(std::ostream& o) { o << color(0); return o; }

ostream& red(ostream& o) {      o << color(31); return o;}
ostream& green(ostream& o) {    o << color(32); return o;}
ostream& brown(ostream& o) {    o << color(33); return o;}
ostream& blue(ostream& o) {     o << color(34); return o;}
ostream& magenta(ostream& o) {  o << color(35); return o;}
ostream& cyan(ostream& o) {     o << color(36); return o;}
ostream& gray(ostream& o) {     o << color(37); return o;}

 /** Bright/Bold **/

ostream& Red(ostream& o) {      o <<color(91); return o;}
ostream& Green(ostream& o) {    o <<color(92); return o;}
ostream& Brown(ostream& o) {    o <<color(93); return o;}
ostream& Blue(ostream& o) {     o <<color(94); return o;}
ostream& Magenta(ostream& o) {  o <<color(95); return o;}
ostream& Cyan(ostream& o) {     o <<color(96); return o;}
ostream& Gray(ostream& o) {     o << color(97); return o;}
