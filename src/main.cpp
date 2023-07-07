#include "gl/glstate.h"
#include "gl/marching_cubes.h"
#include "gl/render.h"
#include "gl/vertex_buffer.h"
#include "gl/window.h"
#include "hf/basis.h"
#include "hf/chem.h"
#include "hf/hf.h"
#include "include/stacktrace.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <eigen3/Eigen/Eigenvalues>
#include <fmt/ranges.h>
#include <fstream>
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_SOURCE[] = {{GL_DEBUG_SOURCE_API, "API"},
                                                                            {GL_DEBUG_SOURCE_WINDOW_SYSTEM, "Window System"},
                                                                            {GL_DEBUG_SOURCE_SHADER_COMPILER, "Shader Compiler"},
                                                                            {GL_DEBUG_SOURCE_THIRD_PARTY, "Third Party"},
                                                                            {GL_DEBUG_SOURCE_APPLICATION, "Application"},
                                                                            {GL_DEBUG_SOURCE_OTHER, "Other"}};

inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_TYPE[] = {{GL_DEBUG_TYPE_ERROR, "Error"},
                                                                          {GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR, "Deprecated Behavior"},
                                                                          {GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR, "Undefined Behavior"},
                                                                          {GL_DEBUG_TYPE_PORTABILITY, "Portability"},
                                                                          {GL_DEBUG_TYPE_PERFORMANCE, "Performance"},
                                                                          {GL_DEBUG_TYPE_MARKER, "Marker"},
                                                                          {GL_DEBUG_TYPE_PUSH_GROUP, "Push Group"},
                                                                          {GL_DEBUG_TYPE_POP_GROUP, "Pop Group"},
                                                                          {GL_DEBUG_TYPE_OTHER, "Other"}};

inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_LEVEL[] = {{GL_DEBUG_SEVERITY_HIGH, "High"},
                                                                           {GL_DEBUG_SEVERITY_MEDIUM, "Medium"},
                                                                           {GL_DEBUG_SEVERITY_LOW, "Low"},
                                                                           {GL_DEBUG_SEVERITY_NOTIFICATION, "Notification"}};

auto string_lookup(GLenum value, const auto& stringList) -> const char*
{
    for (const auto& pair : stringList)
    {
        if (pair.first == value)
        {
            return pair.second;
        }
    }
    return "Unknown";
}

void APIENTRY debugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message,
                                   const void* userParam)
{
    return;
    const char* src = string_lookup(source, std::span(GL_DEBUG_SOURCE));
    const char* ty = string_lookup(type, GL_DEBUG_TYPE);
    const char* level = string_lookup(severity, GL_DEBUG_LEVEL);

    std::cout << "OpenGL Debug Message:" << std::endl;
    std::cout << "  Source: " << src << std::endl;
    std::cout << "  Type: " << ty << std::endl;
    std::cout << "  ID: " << id << std::endl;
    std::cout << "  Severity: " << level << std::endl;
    std::cout << "  Message: " << message << std::endl;

    stacktrace::dump_stacktrace();

    std::cout << std::endl;
}

class mol_orbital : public gl::window
{
    std::vector<hf_isosurface> iso;

    using atoms_vertex_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{
        .enable_normal = true,
        .enable_color = true,
    }>;

    using bound_box_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{.enable_normal = true, .enable_color = true, .primitive = GL_LINES}>;

    atoms_vertex_buffer atoms_vb;
    bound_box_buffer bounding_box_vb;

    hartree_fock_result result;
    gl::render_manager man;

    int index = 0;
    int rot_x = 0;
    int rot_y = 0;
    glm::quat rot = glm::quat(1, 0, 0, 0);
    double zoom = 0;

    int mouse_x = 0;
    int mouse_y = 0;
    bool mouse_init = false;
    bool display_mo = true;
    bool display_atoms = true;
    bool is_left_pressed = false;

    inline static constexpr glm::dvec3 CUBE_CORNER_1 = {-6, -6, -6};
    inline static constexpr glm::dvec3 CUBE_CORNER_2 = {6, 6, 6};

    void next_mo()
    {
        index++;
        index = (index + result.orbitals.size()) % result.orbitals.size();
    }

    void prev_mo()
    {
        index--;
        index = (index + result.orbitals.size()) % result.orbitals.size();
    }

    void regenerate_surfaces()
    {
        auto start = std::chrono::high_resolution_clock::now();
        iso.resize(result.orbitals.size());
        std::cout << "Generating your MOs\n";
        for (size_t i = 0; i < result.orbitals.size(); i++)
        {
            std::cout << "MO: " << i << '\n';
            iso[i].isolevel_hf(result, i, CUBE_CORNER_1, CUBE_CORNER_2, 20);
        }
        auto stop = std::chrono::high_resolution_clock::now();

        std::cout << "Time taken by function: " << (duration_cast<std::chrono::microseconds>(stop - start)).count() << " microseconds" << std::endl;
    }

    void draw_sphere(float r, glm::vec3 center, uint32_t color)
    {
        static constexpr float di = 0.02;
        static constexpr float dj = 0.04;
        float db = di * 2 * M_PI;
        float da = dj * M_PI;

        for (float i = 0; i < 1.0; i += di)
        {
            for (float j = 0; j < 1.0; j += dj)
            {
                float b = i * 2 * M_PI;
                float a = (j - 0.5) * M_PI;

                gl::vec3 normal = {std::cos(a + da / 2) * std::cos(b + db / 2), std::cos(a + da / 2) * std::sin(b + db / 2), std::sin(a + da / 2)};
                gl::vec3 p0(r * std::cos(a) * std::cos(b), r * std::cos(a) * std::sin(b), r * std::sin(a));
                gl::vec3 p1(r * std::cos(a) * std::cos(b + db), r * std::cos(a) * std::sin(b + db), r * std::sin(a));
                gl::vec3 p2(r * std::cos(a + da) * std::cos(b + db), r * std::cos(a + da) * std::sin(b + db), r * std::sin(a + da));
                gl::vec3 p3(r * std::cos(a + da) * std::cos(b), r * std::cos(a + da) * std::sin(b), r * std::sin(a + da));

                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p0 + center).color(color).normal(normal).end());
                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p1 + center).color(color).normal(normal).end());
                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p2 + center).color(color).normal(normal).end());

                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p2 + center).color(color).normal(normal).end());
                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p3 + center).color(color).normal(normal).end());
                atoms_vb.vert(atoms_vertex_buffer::builder().vert(p0 + center).color(color).normal(normal).end());
            }
        }
    }

#define select(n, v) ((EDGES[i] & (1 << (n))) ? (CUBE_CORNER_1.v) : (CUBE_CORNER_2.v))

    void generate_atoms()
    {
        // reset the vb
        atoms_vb = atoms_vertex_buffer();
        for (const auto& i : result.atoms)
        {
            draw_sphere(0.2, i.second, ELEMENT_COLORS[i.first - 1]);
        }
        atoms_vb.bake();

        static constexpr uint8_t EDGES[12] = {0b000100, 0b100110, 0b110010, 0b010000, 0b001101, 0b101111,
                                              0b111011, 0b011001, 0b000001, 0b100101, 0b110111, 0b010011};

        for (uint8_t i = 0; i < 12; i++)
        {
            bounding_box_vb.vert(bound_box_buffer::builder().vert(select(0, x), select(1, y), select(2, z)).color(0xc8c8ff).end());
            bounding_box_vb.vert(bound_box_buffer::builder().vert(select(3, x), select(4, y), select(5, z)).color(0xc8c8ff).end());
        }
        bounding_box_vb.bake();
    }

protected:
    void init() override
    {
        glDebugMessageCallback(debugMessageCallback, nullptr);
        gl::gl_check();
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
        gl::gl_check();
        fov = 5;
        z_near = 0.1;
        z_far = 100;
        bg = {0.1, 0.1, 0.1, 1};
        man.get_shader_manager()
            .add_shader("main", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl")
            .compile_shader("main")
            .use_shader("main");
        regenerate_surfaces();
        generate_atoms();
    }

    void controls(int key, int /*scancode*/, int action, int /*mods*/) override
    {
        if (action == GLFW_PRESS)
        {
            switch (key)
            {
            case GLFW_KEY_LEFT:
                prev_mo();
                break;
            case GLFW_KEY_RIGHT:
                next_mo();
                break;
            case GLFW_KEY_ESCAPE:
                close();
                break;
            case GLFW_KEY_R:
                regenerate_surfaces();
                break;
            case GLFW_KEY_W:
                rot_x++;
                break;
            case GLFW_KEY_A:
                rot_y--;
                break;
            case GLFW_KEY_S:
                rot_x--;
                break;
            case GLFW_KEY_D:
                rot_y++;
                break;
            case GLFW_KEY_Z:
                display_mo = !display_mo;
                break;
            case GLFW_KEY_X:
                display_atoms = !display_atoms;
                break;
            }
        }
    }

    void mause(double xpos, double ypos) override
    {
        if (mouse_init && is_left_pressed)
        {
            auto vec = glm::vec2(xpos - mouse_x, ypos - mouse_y);
            auto axis = glm::normalize(glm::vec3(vec.x, vec.y, 0));
            auto quat = glm::angleAxis(vec.length() * 0.01f, axis);
            rot = glm::normalize(rot * quat);
        }

        mouse_x = xpos;
        mouse_y = ypos;
        mouse_init = true;
    }

    void mause_button(int button, int action, int /*mods*/) override
    {
        if (action == GLFW_PRESS)
        {
            switch (button)
            {
            case GLFW_MOUSE_BUTTON_LEFT:
                is_left_pressed = true;
                break;
            case GLFW_MOUSE_BUTTON_4:
                next_mo();
                break;
            case GLFW_MOUSE_BUTTON_5:
                prev_mo();
                break;
            }
        }
        else
        {
            switch (button)
            {
            case GLFW_MOUSE_BUTTON_LEFT:
                is_left_pressed = false;
                break;
            }
        }
    }

    void mause_scroll(double x_off, double y_off) override { zoom += y_off; }

    void draw_atoms()
    {
        if (!man.get_shader_manager().has_shader("atoms_shader"))
        {
            man.get_shader_manager()
                .add_shader("atoms_shader", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl")
                .compile_shader("atoms_shader");
        }

        man.get_shader_manager()
            .push_shader("atoms_shader")
            .set_vec_float("diffuse_pos", {10, 10, 10})
            .set_vec_float("diffuse_color", {1, 1, 1})
            .set_vec_float("ambient_color", {1, 1, 1})
            .set_float("ambient_strength", 0.2);

        man.render_vertex_buffer(atoms_vb);

        man.get_shader_manager()
            .push_shader("atoms_shader")
            .set_vec_float("diffuse_pos", {10, 10, 10})
            .set_vec_float("diffuse_color", {0, 0, 0})
            .set_vec_float("ambient_color", {1, 1, 1})
            .set_float("ambient_strength", 1);

        man.render_vertex_buffer(bounding_box_vb);

        man.get_shader_manager().pop_shader();
    }

    void render(const gl::mat4& proj) override
    {
        man.view().push().mult(glm::lookAt(glm::vec3(10 * std::exp(zoom * 0.693), 0, 0), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0)));
        man.projection().push().mult(proj);
        man.model().push().rotate(-rot_x * 0.1F, 0, 0, 1).rotate(rot_y * 0.1F, 0, 1, 0);

        gl::gl_check();

        if (display_mo)
        {
            iso[index].draw(man);
        }

        if (display_atoms)
        {
            draw_atoms();
        }

        man.model().pop();
        man.projection().pop();
        man.view().pop();
    }

public:
    mol_orbital(int a, int b, const char* c, hartree_fock_result result) : gl::window(a, b, c), result(std::move(result)) {}
};

auto main(int argc, char* argv[]) -> int
{
    try
    {
        std::string filename;

        if (argc > 1)
        {
            filename = argv[argc - 1];
        }

        basis_manager::get_instance().register_basis("sto-2g");
        basis_manager::get_instance().register_basis("sto-3g");
        basis_manager::get_instance().register_basis("sto-4g");
        basis_manager::get_instance().register_basis("sto-5g");
        basis_manager::get_instance().register_basis("sto-6g");

        auto mol = molecule::read_from_file(filename);
        auto res = solve(mol);

        std::cout << fmt::format("iter={}, energies={}\n", res.iterations,
                                 res.mo_energies | std::views::transform([](double x) { return x * 27.2114; }));
        mol_orbital w(120, 120, "Molecular Orbitals", std::move(res));

        w.run();
        return 0;
    }
    catch (gl::shader_error& err)
    {
        std::cerr << fmt::format("unhandled shader exception was received: {}\n", err.what());
        std::cerr << fmt::format("logs: \n{}\n", err.get_log());
        exit(-1);
    }
    catch (std::exception& err)
    {
        std::cerr << fmt::format("unhandled exception was received: {}\n", err.what());
        exit(-1);
    }
    if (argc == 1)
    {
        std::cout << "Usage: " << argv[0] << " <input file>" << std::endl;
        std::cout << "Please provide a valid input file" << std::endl;
        return 0;
    }
}

