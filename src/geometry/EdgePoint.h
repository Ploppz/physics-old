#include <iostream>
#include <glm/glm.hpp>

// External
class Polygon;


struct EdgePoint {
    EdgePoint() {};
    EdgePoint(Polygon* parent) : parent(parent){};
    EdgePoint(int index, float alpha, Polygon* parent);
    void set(int index, float alpha) {this->index = index; this->alpha = alpha; }
    void set(int index, float alpha, Polygon* parent) {this->index = index; this->alpha = alpha; this->parent = parent; }

    glm::vec2 point();
    glm::vec2 point_t(); //transformed point

    int index;
    float alpha;
    Polygon* parent;

    const bool operator== (const EdgePoint& other);
};
std::ostream& operator<< (std::ostream& o, EdgePoint p);
