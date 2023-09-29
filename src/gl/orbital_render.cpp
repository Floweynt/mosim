#include "gl/orbital_render.h"
#include "gl/glstate.h"
#include "gl/shader.h"
#include "hf/hf.h"
#include <fmt/ranges.h>

struct GTO
{
    float coeff_norm;
    float alpha;
    float x;
    float y;
    float z;
    uint32_t l;
    uint32_t m;
    uint32_t n;
};

inline static constexpr auto BUFFER_SIZE = 1024 * 1024 * 4;

void hf_isosurface::isolevel_hf(const hartree_fock_result& result, uint32_t which_mo, glm::vec3 start, glm::vec3 end, size_t detail, float isolevel)
{
    static gl::compute_shader shader("assets/shaders/render_cloud_compute.glsl");
    if (!shader.get_id())
    {
        shader.compile();
    }

    shader.set_float("isolevel", isolevel);
    shader.set_vec_float("start", start);
    shader.set_vec_float("cube_step", (end - start) / float(4 * detail));
    shader.set_uint("orbital_count", result.orbitals.size());
    shader.set_uint("gto_per_orbital", result.orbitals[0].get_data().orbs().size());

    // set output buffer
    vbo.init_empty(BUFFER_SIZE, GL_STATIC_DRAW);
    vbo.write(0, uint32_t(0));

    // setup coefficient buffer
    gl::buffer_object mo_coeff_buffer;
    std::vector<float> coeff(result.coefficients.cols());
    for (size_t i = 0; i < result.coefficients.cols(); i++)
    {
        coeff[i] = (float)result.coefficients(i, which_mo);
    }
    mo_coeff_buffer.init(std::span(coeff), GL_STATIC_DRAW);

    // setup GTO buffer
    gl::buffer_object gtos_buffer;
    std::vector<GTO> gtos;
    for (const auto& cgf : result.orbitals)
    {
        for (const auto& gto : cgf.orbs())
        {
            gtos.emplace_back(gto.factor(), gto.get_alpha(), cgf.get_pos().x, cgf.get_pos().y, cgf.get_pos().z, gto.get_l(), gto.get_m(),
                              gto.get_n());
        }
    }

    gtos_buffer.init(std::span(gtos), GL_STATIC_DRAW);

    shader.use();

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, mo_coeff_buffer.get_id());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, gtos_buffer.get_id());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, vbo.get_id());
    glDispatchCompute(detail, detail, detail);
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, 0);

    vert_count = vbo.read<uint32_t>(0) / 2;
    vao.bind_vbo(vbo, 0, 4 * sizeof(float), 8 * sizeof(float));
    vao.format(0, 4, GL_FLOAT, false, 0);
    vao.enable_attribute(0);
    vao.set_attribute_buffer_bind(0, 0);
    vao.format(1, 4, GL_FLOAT, false, 4 * sizeof(float));
    vao.enable_attribute(1);
    vao.set_attribute_buffer_bind(1, 0);
}

void hf_isosurface::render(gl::render_manager& render, glm::vec3 light_pos)
{
    if (!render.get_shader_manager().has_shader("__hf_isosurface"))
    {
        render.get_shader_manager()
            .add_shader("__hf_isosurface", "assets/shaders/orbital_vertex.glsl", "assets/shaders/orbital_fragment.glsl",
                        "assets/shaders/orbital_geom.glsl")
            .compile_shader("__hf_isosurface");
    }

    render.get_shader_manager()
        .push_shader("__hf_isosurface")
        .set_vec_float("diffuse_pos", light_pos)
        .set_vec_float("diffuse_color", {1, 1, 1})
        .set_vec_float("ambient_color", {1, 1, 1})
        .set_float("ambient_strength", 0.5);

    render.render_generic([this]() { vao.draw(GL_POINTS, 0, vert_count); });
    render.get_shader_manager().pop_shader();
}

