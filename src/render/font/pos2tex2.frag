#version 150

in vec2 f_uv;
in vec3 f_color;
uniform sampler2D tex;
out vec4 outColor;

float limit = 0;

float sampl(in vec2 uv, float width)
{
	float sdf = texture(tex, uv).r;
	return smoothstep(limit - width, limit + width, sdf);
}

void main()
{
	float data = texture(tex, f_uv).r;
	float width = fwidth(data);

	/* float alpha = sampl(f_uv, width); */
	float alpha = texture(tex, f_uv).r;

	/* outColor = vec4(0, 0, 0, mix(1, 0, alpha)); */
	outColor = vec4(vec3(f_color), alpha);
}
