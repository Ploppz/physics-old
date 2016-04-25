#include "tmp.h"
#include "glutils.h"

#include <glm/glm.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>


// Uses Polygon p as reference coordinate system
unsigned int addToBuffer(LineSegment e, float *buffer, int offset)
{
	buffer[offset] = e.first.x; ++ offset;
	buffer[offset] = e.first.y; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 0; ++ offset;

	buffer[offset] = e.second.x; ++ offset;
	buffer[offset] = e.second.y; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 1; ++ offset;
	buffer[offset] = 0; ++ offset;
	return offset;
}

void addToBuffer(LineSegment e, std::vector<float> &buffer)
{
	buffer.push_back(e.first.x);
	buffer.push_back(e.first.y);
	buffer.push_back(1);
	buffer.push_back(1);
	buffer.push_back(0);

	buffer.push_back(e.second.x);
	buffer.push_back(e.second.y);
	buffer.push_back(1);
	buffer.push_back(1);
	buffer.push_back(0);
}
void addToBuffer(LineSegment e, BufferWriter<float> & buffer, float r, float g, float b)
{
	buffer.write(e.first.x, e.first.y);
    buffer.write(r, g, b);

	buffer.write(e.second.x, e.second.y);
    buffer.write(r, g, b);
}

float randFloat() {
    return rand() / float(RAND_MAX);
}

glm::vec3 randomColor()
{
    return glm::vec3(randFloat(), randFloat(), randFloat()) * 0.5f;
}

bool file_exists(const std::string &name)
{
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}   
}
void assert_not_equal(glm::vec2 a, glm::vec2 b)
{
    assert(a.x != b.x && a.y != b.y);
}
