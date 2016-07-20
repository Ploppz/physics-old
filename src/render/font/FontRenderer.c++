#include "FontRenderer.h"
#include "shaders.h"
#include "glutils.h"
#include "tmp.h"

/* glm */
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

/* std */
#include <string>
#include <iostream>
#include <algorithm>

const int GL_VERTEX_STRIDE = 7;

FontRenderer::FontRenderer(int size, const std::string font_file_name)
	: size(size), color(1)
{
	std::string s;
	if (file_exists("font_atlas")) {
		texture.useExistingAtlas("font_atlas", "metadata_atlas");
	} else {
		texture.generateAtlas(font_file_name);
		texture.writeAtlasToFile("font_atlas", "metadata_atlas");
	}
    setup_GL();
}
FontRenderer::FontRenderer()
    : FontRenderer(1, "fonts/peep-07x14.bdf")
{ }


void FontRenderer::add_text(const std::string& text, float penx, float peny, bool kerning) {
    add_text(text, penx, peny, 1, kerning);
}
void FontRenderer::add_text(const std::string& text, float penx, float peny, float scale, bool kerning)
{
	// x & y are pen coordinates
	Glyph g;
	unsigned char prevChar = 0;
	for (auto it = text.begin(); it != text.end(); it ++)
	{
		g = texture.getGlyph(*it);
		if (prevChar && kerning) {
			penx += g.kerning[prevChar] * size;
			// std::cout << "Kerning for " << *it << ": " << g.kerning[prevChar] * size << std::endl;
		}
		// x and y, offset
		float x = penx + g.xoffset * size;
		float y = peny - g.yoffset * size;

		// std::cout << "offset for " << *it << ": (" << g.xoffset << ", " << g.yoffset << ")" << std::endl;
		
		// Actual width and height
		float w = g.width * size * scale;
		float h = g.height * size * scale;
		// Normalized texture coordinates
		float u = static_cast<float>(g.u) / texture.getAtlasWidth();
		float v = static_cast<float>(g.v) / texture.getAtlasHeight();
		// Normalized texture width and height
		float tw = static_cast<float>(g.twidth) / texture.getAtlasWidth();
		float th = static_cast<float>(g.theight) / texture.getAtlasHeight();

		active_buffer->push_back(x);
		active_buffer->push_back(y);
		active_buffer->push_back(u);
		active_buffer->push_back(v);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		active_buffer->push_back(x + w);
		active_buffer->push_back(y);
		active_buffer->push_back(u + tw);
		active_buffer->push_back(v);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		active_buffer->push_back(x + w);
		active_buffer->push_back(y - h);
		active_buffer->push_back(u + tw);
		active_buffer->push_back(v + th);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		active_buffer->push_back(x);
		active_buffer->push_back(y);
		active_buffer->push_back(u);
		active_buffer->push_back(v);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		active_buffer->push_back(x);
		active_buffer->push_back(y - h);
		active_buffer->push_back(u);
		active_buffer->push_back(v + th);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		active_buffer->push_back(x + w);
		active_buffer->push_back(y - h);
		active_buffer->push_back(u + tw);
		active_buffer->push_back(v + th);
        active_buffer->push_back(color.r); active_buffer->push_back(color.g); active_buffer->push_back(color.b);

		//TODO: Kerning
		penx += g.xadvance * size * scale;
		// std::cout << "Advance for " << *it << ": " << g.xadvance * size << std::endl;

		prevChar = *it;
	}
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

void FontRenderer::setup_GL()
{
	GLuint vertexShader, fragmentShader; // unused
	shader_program  = createShaderProgram(shaders::render_font_pos2tex2_v, shaders::render_font_pos2tex2_f, vertexShader, fragmentShader);
	glUseProgram(shader_program);


    // Set initial size of buffer
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, 8000 * sizeof(float), NULL, GL_STREAM_DRAW);

	vao = createVertexArrayObject();
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	setFormat("position 2f texcoor 2f color 3f", shader_program);

	// Textures
	glGenTextures(1, &texture_ref);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texture_ref);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, texture.getAtlasWidth(), texture.getAtlasHeight()
			, 0, GL_RED, GL_FLOAT, texture.atlas.data());
	glUniform1i(glGetUniformLocation(shader_program, "tex"), 0);
	// Mipmap
	glGenerateMipmap(GL_TEXTURE_2D);
	//
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}
void FontRenderer::render(float center_x, float center_y, int width, int height, float zoom)
{
    if (active_buffer->size() == 0) return;
	glUseProgram(shader_program);
	glBindVertexArray(vao);
	glBindTexture(GL_TEXTURE_2D, texture_ref);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	// Reupload
	float *buffer_ptr = static_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, active_buffer->size() * sizeof(float), GL_MAP_WRITE_BIT));
    std::copy_n(active_buffer->data(), active_buffer->size(), buffer_ptr);

    glUnmapBuffer(GL_ARRAY_BUFFER);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnablei(GL_BLEND, 0);
	// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	// Upload matrices
	glm::mat4 proj = ortho2D(width * zoom, height * zoom, 0, 1);
	glm::mat4 view = viewMatrix2D(center_x, center_y, 1, 1);
	glm::mat4 model = glm::mat4();
	GLuint projUni = glGetUniformLocation(shader_program, "proj");
	GLuint viewUni = glGetUniformLocation(shader_program, "view");
	GLuint modelUni = glGetUniformLocation(shader_program, "model");
	glUniformMatrix4fv(projUni, 1, GL_FALSE, glm::value_ptr(proj));
	glUniformMatrix4fv(viewUni, 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(modelUni, 1, GL_FALSE, glm::value_ptr(model));

	// Draw
	glDrawArrays(GL_TRIANGLES, 0, active_buffer->size() / GL_VERTEX_STRIDE);
}
void FontRenderer::set_color(float r, float g, float b)
{
    color = glm::vec3(r,g,b);
}
