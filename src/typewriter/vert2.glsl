#version 150

in vec2 position;
in vec3 color;

out vec4 Color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;
uniform float scale;

void main() {
	Color = vec4(color, 0.5);
	
	gl_Position = vec4(position, 0, 1);
}
