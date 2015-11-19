#include "Polygon.h"
#include "BodySystem.h"
#include "glutils.h"

#include <vector>
#include <glm/glm.hpp>

// Temporary place for code


unsigned int addToBuffer(LineSegment e, Body &b, float *buffer, int offset);
void addToBuffer(LineSegment e, Body &b, std::vector<float> &buffer);
void addToBuffer(LineSegment e, Body &b, BufferWriter<float> &buffer);

glm::vec3 randomColor();

template <typename T>
class Graph
{
public:
    class Node;
    std::vector<Node> nodes;

    class Node
    {
        T value;
        std::vector<int> neighbors; // indices into graph
    };
private:
};
