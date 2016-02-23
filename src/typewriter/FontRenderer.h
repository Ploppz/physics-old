#include "FontTexture.h"
#include <string>
#include <GL/glew.h>
#include <glm/glm.hpp>

struct FontVertex
{
	//TODO: Add z
	FontVertex(float x, float y, float z, float u, float v) {
		this->x = x; this->y = y; this->u = u; this->v = v;
	}
	// position 2f texcoor 2f
	float x, y, u, v;
};

class FontRenderer
{
private:
	int size;
	// Structure: position 2f texcoor 2f
	std::vector<float> buffer;
	FontTexture &texture;
    glm::vec3 color;


	// OpenGL
	GLuint shaderProgram;
	GLuint textureRef;
	GLuint bufferRef;
	GLuint VAO;

public:
	FontRenderer(int size, FontTexture &texture);
	
    void setColor(int r, int g, int b);
	void addText(const std::string text, float x, float y, bool kerning);
	void clearBuffer();
	void setup();
	void render(float width, float height);
};
