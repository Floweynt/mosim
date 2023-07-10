#version 330 core

layout(location = 0) in vec4 in_pos;
layout(location = 1) in vec4 in_norm;
// layout (location = 2) in vec2 in_uv;
// layout (location = 3) in uint in_color;
// layout (location = 3) in vec4 in_color;

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;

out vec3 vertex_normal;
out vec4 vertex_pos;

void main()
{
    gl_Position = (projection * view * model) * vec4(in_pos.xyz, 1.0);
    vertex_normal = mat3(transpose(inverse(model))) * vec3(in_norm);
    
    vec4 in_pos_real = vec4(in_pos.xyz, 1);
    float in_pos_phase = in_pos.w;
    vertex_pos = vec4((model * in_pos_real).xyz, in_pos_phase);
}
