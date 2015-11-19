#version 150

in vec2 position;
in vec2 texcoor;

out vec2 Texcoor;
out vec2 Pos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;
uniform float scale;

void main() {
	Texcoor.x = 0.5 + (texcoor.x - 0.5) / scale;
	Texcoor.y = 0.5 + (texcoor.y - 0.5) / scale;
	Pos = position;
	gl_Position = vec4(position, 0, 1);
}
