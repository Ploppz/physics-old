#version 150

uniform mat4 proj, view, model;
in vec2 position;

void main()
{
	gl_Position = proj * view * model * vec4(position, 0, 1);
}
