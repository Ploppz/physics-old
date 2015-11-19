#include "tmp.h"
#include "glutils.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>


// Uses Polygon p as reference coordinate system
unsigned int addToBuffer(LineSegment e, Body &b, float *buffer, int offset)
{
	buffer[offset] = e.first.x + b.position().x; ++ offset;
	buffer[offset] = e.first.y + b.position().y; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 0; ++ offset;

	buffer[offset] = e.second.x + b.position().x; ++ offset;
	buffer[offset] = e.second.y + b.position().y; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 0; ++ offset;
	return offset;
}

void addToBuffer(LineSegment e, Body &b, std::vector<float> &buffer)
{
	buffer.push_back(e.first.x + b.position().x);
	buffer.push_back(e.first.y + b.position().y);
	buffer.push_back(1);
	buffer.push_back(1);
	buffer.push_back(0);

	buffer.push_back(e.second.x + b.position().x);
	buffer.push_back(e.second.y + b.position().y);
	buffer.push_back(1);
	buffer.push_back(1);
	buffer.push_back(0);
}
void addToBuffer(LineSegment e, Body &b, BufferWriter<float> & buffer)
{
	buffer.write(e.first.x + b.position().x, e.first.y + b.position().y);
    buffer.write(1, 1, 0);

	buffer.write(e.second.x + b.position().x, e.second.y + b.position().y);
    buffer.write(1, 1, 0);
}

float randFloat() {
    return rand() / float(RAND_MAX) / 2;
}

glm::vec3 randomColor()
{
    return glm::vec3(randFloat(), randFloat(), randFloat());
}

