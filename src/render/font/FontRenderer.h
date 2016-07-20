#pragma once
#include "FontTexture.h"
#include "render/Renderer.h"

#include <string>
#include <GL/glew.h>
#include <glm/glm.hpp>


// Buffer format: position 2f texcoor 2f


class FontRenderer : public Renderer
{

 public:
    FontRenderer();
	FontRenderer(int size, const std::string font_file_name);
	void render(float center_x, float center_y, int width, int height, float zoom);
    void set_color(float r, float g, float b);
    // Draw //
	void add_text(const std::string& text, float x, float y, bool kerning);
    void add_text(const std::string& text, float penx, float peny, float scale, bool kerning);
 private:
	void setup_GL();
 private:
	int size;
	FontTexture texture;
    glm::vec3 color;

	// GL //
	GLuint shader_program;
	GLuint texture_ref;
	GLuint vbo;
	GLuint vao;
};
