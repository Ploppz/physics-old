#version 150

uniform mat4 proj, view, model;
in vec2 position;
in vec3 color;

out vec3 Color;

void main()
{
	Color = color;
	gl_Position = proj * view * model * vec4(position, 0, 1);
}
