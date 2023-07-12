#version 330 core
layout (location = 0) in vec3 vertex;
layout (location = 3) in uint in_color;

out vec4 vert_color;

uniform mat4 projection;

void main()
{
    float a = float((in_color & 0xff000000u) >> 24) / 255.f;
    float r = float((in_color & 0x00ff0000u) >> 16) / 255.f;
    float g = float((in_color & 0x0000ff00u) >> 8 ) / 255.f;
    float b = float((in_color & 0x000000ffu)      ) / 255.f;

    vert_color = vec4(r, g, b, 1);
    gl_Position = projection * vec4(vertex, 1.0);
}


