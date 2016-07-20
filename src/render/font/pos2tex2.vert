#version 150

in vec2 position;
in vec2 texcoor;
in vec3 color;

out vec2 f_uv;
out vec3 f_color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main()
{
	f_uv = texcoor;
    f_color = color;
	gl_Position = proj * view * model * vec4(position, 0, 1);
}
