#version 330 core

layout (location = 0) in vec3 in_pos;
layout (location = 1) in vec3 in_norm;
// layout (location = 2) in vec2 in_uv;
// layout (location = 3) in uint in_color;
layout (location = 3) in uint in_color;

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;

out vec4 out_color;
out vec3 out_normal;
out vec3 out_pos;

void main()
{
    gl_Position = (projection * view * model) * vec4(in_pos.xyz, 1.0);
    float a = float((in_color & 0xff000000u) >> 24) / 255.f;
    float r = float((in_color & 0x00ff0000u) >> 16) / 255.f;
    float g = float((in_color & 0x0000ff00u) >> 8 ) / 255.f;
    float b = float((in_color & 0x000000ffu)      ) / 255.f;

    out_color = vec4(r, g, b, 1);
    out_normal = mat3(transpose(inverse(model))) * vec3(in_norm);
    out_pos = vec3(model * vec4(in_pos, 1));
}
