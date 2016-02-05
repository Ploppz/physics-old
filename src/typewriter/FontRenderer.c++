#include "FontRenderer.h"
#include "../shaders.h"
#include "../glutils.h"

/* glm */
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

/* std */
#include <string>
#include <iostream>
#include <algorithm>

FontRenderer::FontRenderer(int size, FontTexture &texture)
	: size(size), texture(texture)
{}

void FontRenderer::addText(std::string text, float penx, float peny, bool kerning)
{
	// x & y are pen coordinates
	Glyph g;
	unsigned char prevChar = 0;
	float x, y, u, v, w, h, tw, th;
	for (auto it = text.begin(); it != text.end(); it ++)
	{
		g = texture.getGlyph(*it);
		if (prevChar && kerning) {
			penx += g.kerning[prevChar] * size;
			// std::cout << "Kerning for " << *it << ": " << g.kerning[prevChar] * size << std::endl;
		}
		// x and y, offset
		x = penx + g.xoffset * size;
		y = peny - g.yoffset * size;

		// std::cout << "offset for " << *it << ": (" << g.xoffset << ", " << g.yoffset << ")" << std::endl;
		
		// Actual width and height
		w = g.width * size;
		h = g.height * size;
		// Normalized texture coordinates
		u = static_cast<float>(g.u) / texture.getAtlasWidth();
		v = static_cast<float>(g.v) / texture.getAtlasHeight();
		// Normalized texture width and height
		tw = static_cast<float>(g.twidth) / texture.getAtlasWidth();
		th = static_cast<float>(g.theight) / texture.getAtlasHeight();

		buffer.push_back(x);
		buffer.push_back(y);
		buffer.push_back(u);
		buffer.push_back(v);

		buffer.push_back(x + w);
		buffer.push_back(y);
		buffer.push_back(u + tw);
		buffer.push_back(v);

		buffer.push_back(x + w);
		buffer.push_back(y - h);
		buffer.push_back(u + tw);
		buffer.push_back(v + th);

		buffer.push_back(x);
		buffer.push_back(y);
		buffer.push_back(u);
		buffer.push_back(v);

		buffer.push_back(x);
		buffer.push_back(y - h);
		buffer.push_back(u);
		buffer.push_back(v + th);

		buffer.push_back(x + w);
		buffer.push_back(y - h);
		buffer.push_back(u + tw);
		buffer.push_back(v + th);

		//TODO: Kerning
		penx += g.xadvance * size;
		// std::cout << "Advance for " << *it << ": " << g.xadvance * size << std::endl;

		prevChar = *it;
	}
}
void FontRenderer::clearBuffer()
{
	buffer.clear();
}

GLfloat rectangle1[] = {
	// x, y, tex_x, tex_y
	-1.0f, -1.0f, 0.0f, 1.0f,
	1.0f, -1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f, 0.0f,

	-1.0f, -1.0f,	0.0f, 1.0f, 
	1.0f, 1.0f, 1.0f, 0.0f,
	-1.0f, 1.0f, 0.0f, 0.0f
};

void FontRenderer::setup()
{
	GLuint vertexShader, fragmentShader; // unused
	shaderProgram  = createShaderProgram(shaders::typewriter_pos2tex2_v, shaders::typewriter_pos2tex2_f, vertexShader, fragmentShader);
	glUseProgram(shaderProgram);


    // Set initial size of buffer
	glGenBuffers(1, &bufferRef);
	glBindBuffer(GL_ARRAY_BUFFER, bufferRef);
	glBufferData(GL_ARRAY_BUFFER, 8000 * sizeof(float), NULL, GL_STREAM_DRAW);

	VAO = createVertexArrayObject();
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, bufferRef);
	setFormat("position 2f texcoor 2f", shaderProgram);

	// Textures
	glGenTextures(1, &textureRef);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, textureRef);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, texture.getAtlasWidth(), texture.getAtlasHeight()
			, 0, GL_RED, GL_FLOAT, texture.atlas.data());
	glUniform1i(glGetUniformLocation(shaderProgram, "tex"), 0);
	// Mipmap
	glGenerateMipmap(GL_TEXTURE_2D);
	//
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}
void FontRenderer::render(float width, float height)
{
    if (buffer.size() == 0) return;
	glUseProgram(shaderProgram);
	glBindVertexArray(VAO);
	glBindTexture(GL_TEXTURE_2D, textureRef);
	glBindBuffer(GL_ARRAY_BUFFER, bufferRef);

	// Reupload
	float *buffer_ptr = static_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, buffer.size() * sizeof(float), GL_MAP_WRITE_BIT));
    std::copy_n(buffer.data(), buffer.size(), buffer_ptr);

    glUnmapBuffer(GL_ARRAY_BUFFER);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnablei(GL_BLEND, 0);
	// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	// Upload matrices
	glm::mat4 proj = ortho2D(width, height, 0, 1);
	glm::mat4 view = viewMatrix2D(0, 0, 1, 1);
	glm::mat4 model = glm::mat4();
	GLuint projUni = glGetUniformLocation(shaderProgram, "proj");
	GLuint viewUni = glGetUniformLocation(shaderProgram, "view");
	GLuint modelUni = glGetUniformLocation(shaderProgram, "model");
	glUniformMatrix4fv(projUni, 1, GL_FALSE, glm::value_ptr(proj));
	glUniformMatrix4fv(viewUni, 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(modelUni, 1, GL_FALSE, glm::value_ptr(model));

	// Draw
	glDrawArrays(GL_TRIANGLES, 0, buffer.size() / 4);
}
