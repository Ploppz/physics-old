#pragma once
#include "glutils.h"

#include <vector>
#include <glm/glm.hpp>

// Temporary place for code
// typedef std::pair<glm::vec2, glm::vec2> LineSegment;


// unsigned int addToBuffer(LineSegment e, float *buffer, int offset);
// void addToBuffer(LineSegment e, std::vector<float> &buffer);
// void addToBuffer(LineSegment e, BufferWriter<float> &buffer, float, float, float);

float randFloat();
glm::vec3 randomColor();
bool file_exists(const std::string &name);

void assert_not_equal(glm::vec2, glm::vec2);
