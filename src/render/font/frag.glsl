#version 150

in vec2 Texcoor;
in vec2 Pos;

out vec4 outColor;

uniform sampler2D tex;
uniform float limit;

void main() {

	vec4 data = texture(tex, Texcoor);
	vec4 no_color = vec4(0, 0, 0, 1);
	vec4 yes_color = vec4(1, 1, 1, 1);
	/* float mixAmount = step(limit, data.r);*/
	float delta = fwidth(data.r);
	float mixAmount = smoothstep(limit - delta, limit, data.r);
	outColor = vec4(vec3(mix(no_color, yes_color, mixAmount)), 1);
}
