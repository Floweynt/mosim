#version 430
#define LOCAL_SIZE 4
#extension GL_EXT_gpu_shader4 : enable
#pragma optionNV(unroll all)
// marching cubes stuff
// in
uniform float isolevel;
uniform vec3 start;
uniform vec3 cube_step;
uniform isampler2D TRIANGLE_LOOKUP;
// out
layout(std430, binding = 2) restrict buffer out_vert_t
{
    uint count;
    vec4 vert_list[];
};

// compute shader stuff
layout(local_size_x = LOCAL_SIZE, local_size_y = LOCAL_SIZE, local_size_z = LOCAL_SIZE) in;

// lookup tables
const uvec3 VERTICES[8] =
    uvec3[](uvec3(0, 0, 0), uvec3(1, 0, 0), uvec3(1, 1, 0), uvec3(0, 1, 0), uvec3(0, 0, 1), uvec3(1, 0, 1), uvec3(1, 1, 1), uvec3(0, 1, 1));
const uvec2 EDGES[12] = uvec2[](uvec2(0, 1), uvec2(1, 2), uvec2(2, 3), uvec2(3, 0), uvec2(4, 5), uvec2(5, 6), uvec2(6, 7), uvec2(7, 4), uvec2(0, 4),
                                uvec2(1, 5), uvec2(2, 6), uvec2(3, 7));

// MO specific stuff
struct GTO
{
    float coeff_norm;
    float alpha;
    float x;
    float y;
    float z;
    uint l;
    uint m;
    uint n;
};

uniform uint orbital_count;
uniform uint gto_per_orbital;
layout(std430, binding = 0) buffer mo_coeff_t
{
    float coeff[]; // this will be of length orbital_count
};
layout(std430, binding = 1) buffer gtos_t
{
    GTO gtos[]; // this will be of length orbital_count * gto_per_orbital
};

float fast_pow(float base, uint exp)
{
    float result = 1;
    while (true)
    {
        if ((exp & 1) != 0)
            result *= base;
        exp >>= 1;
        if (exp == 0)
            break;
        base *= base;
    }

    return result;
}

float phi(vec3 r, uint idx)
{
    float sum = 0;

    uint root_index = idx * gto_per_orbital;

    vec3 R = r - vec3(gtos[root_index].x, gtos[root_index].y, gtos[root_index].z);
    float square_mag = dot(R, R);

    for (uint i = 0; i < gto_per_orbital; i++)
    {
        uint gto_idx = root_index + i;
        sum += gtos[gto_idx].coeff_norm * exp(-gtos[gto_idx].alpha * square_mag);
    }

    return sum * fast_pow(R.x, gtos[root_index].l) * fast_pow(R.y, gtos[root_index].m) * fast_pow(R.z, gtos[root_index].n);
}

float mo_phi(vec3 r)
{
    float sum = 0;
    for (uint i = 0; i < orbital_count; i++)
        if (abs(coeff[i]) > 0.01)
            sum += phi(r, i) * coeff[i];
    return sum;
}

#define EPS 0.0001
const vec3 GRAD_DELTA_X = vec3(EPS, 0, 0);
const vec3 GRAD_DELTA_Y = vec3(0, EPS, 0);
const vec3 GRAD_DELTA_Z = vec3(0, 0, EPS);

vec3 mo_grad(vec3 r)
{
    float val = mo_phi(r);
    return vec3((mo_phi(r + GRAD_DELTA_X) - val) / EPS, (mo_phi(r + GRAD_DELTA_Y) - val) / EPS, (mo_phi(r + GRAD_DELTA_Z) - val) / EPS);
}

vec3 cube_pos(uint i) { return start + vec3(gl_GlobalInvocationID + VERTICES[i]) * cube_step; }

float cube_val(uint i) { return mo_phi(cube_pos(i)); }

vec3 vertex_interp(vec3 p0, vec3 p1, float d0, float d1)
{
    float diff = d1 - d0;
    return (p1 - p0) * (isolevel - d0) / diff + p0;
}

int lookup_vert(uint i, uint j) { return texelFetch2D(TRIANGLE_LOOKUP, ivec2(int(j), int(i)), 0).a; }

void emit_array(in vec4 out_buf[32], int number_of_verts)
{
    uint idx = atomicAdd(count, number_of_verts);
    for (int i = 0; i < number_of_verts; i++)
        vert_list[idx + i] = out_buf[i];
}

void march_the_cubes(in float cube_vals[8], float phase)
{
    uint cube_index = 0;
    for (uint i = 0; i < 8; i++)
        cube_index |= uint(cube_vals[i] * phase < isolevel) * (1 << i);

    if (cube_index == 0 || cube_index == 255)
        return;

    vec3 vert_list[12];
    for (uint i = 0; i < 12; i++)
    {
        vert_list[i] = vertex_interp(cube_pos(EDGES[i].x), cube_pos(EDGES[i].y), cube_vals[EDGES[i].x] * phase, cube_vals[EDGES[i].y] * phase);
    }

    vec4 out_buf[32];

    int i = 0;
    for (; lookup_vert(cube_index, i * 3) != -1; i++)
    {
        vec3 norm_approx = -phase * normalize(mo_grad(vec3(vert_list[lookup_vert(cube_index, i * 3 + 0)])));
        // write the values
        out_buf[i * 6 + 0] = vec4(vert_list[lookup_vert(cube_index, i * 3 + 0)], phase);
        out_buf[i * 6 + 1] = vec4(norm_approx, 1);
        out_buf[i * 6 + 2] = vec4(vert_list[lookup_vert(cube_index, i * 3 + 1)], phase);
        out_buf[i * 6 + 3] = vec4(norm_approx, 1);
        out_buf[i * 6 + 4] = vec4(vert_list[lookup_vert(cube_index, i * 3 + 2)], phase);
        out_buf[i * 6 + 5] = vec4(norm_approx, 1);
    }

    emit_array(out_buf, i * 6);
}

void main()
{
    // if (abs(mo_phi(cube_pos(0))) < isolevel / 10)
    //     return;

    float cube_vals[8];
    for (uint i = 0; i < 8; i++)
        cube_vals[i] = cube_val(i);

    march_the_cubes(cube_vals, 1);
    march_the_cubes(cube_vals, -1);
}

