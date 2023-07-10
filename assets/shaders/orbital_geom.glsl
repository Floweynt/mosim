#version 330

layout(points) in;
layout(points, max_vertices = 1) out;

in vec3 vertex_normal[];
in vec4 vertex_pos[];

out vec4 geom_color;
out vec3 geom_normal;
out vec3 geom_pos;

void main(void)
{
    int i;
    vec4 color = vertex_pos[0].w > 0 ? vec4(0, 1, 0, 1) : vec4(1, 0, 0, 1);
    geom_color = color;
    geom_normal = vertex_normal[0];
    geom_pos = vec3(vertex_pos[0]);
    gl_Position = gl_in[i].gl_Position;
    gl_PointSize = 3;
    EmitVertex();
    EndPrimitive();
}

