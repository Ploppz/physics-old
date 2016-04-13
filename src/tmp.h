#include "BodySystem.h"
#include "glutils.h"

#include <vector>
#include <glm/glm.hpp>

// Temporary place for code


unsigned int addToBuffer(LineSegment e, float *buffer, int offset);
void addToBuffer(LineSegment e, std::vector<float> &buffer);
void addToBuffer(LineSegment e, BufferWriter<float> &buffer, float, float, float);

float randFloat();
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
