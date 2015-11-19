#version 150

in vec2 position;
in vec2 texcoor;

out vec2 f_uv;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main()
{
	f_uv = texcoor;
	gl_Position = proj * view * model * vec4(position, 0, 1);
}
