#include "build_config.h"
#include "gl/font.h"
#include "gl/glstate.h"
#include "gl/orbital_render.h"
#include "gl/render.h"
#include "gl/vertex_buffer.h"
#include "gl/window.h"
#include "hf/basis.h"
#include "hf/chem.h"
#include "hf/hf.h"
#include "include/stacktrace.h"
#include "platform_dep.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fmt/ranges.h>
#include <fstream>
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

class fps_calculator
{
    double reported_fps = 0;
    double reported_min_fps = 0;
    double reported_max_fps = 0;

    int64_t prev_report_tick = 0;
    int64_t prev_tick = 0;
    int64_t frame_count = 0;
    int64_t max_ms = INT_MIN;
    int64_t min_ms = INT_MAX;

public:
    void report_frame()
    {
        int64_t current_tick = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

        // probably not the unix epoch right now...
        if (prev_tick == 0)
        {
            prev_report_tick = prev_tick = current_tick;
            return;
        }

        // update everything
        frame_count++;
        max_ms = std::max(max_ms, current_tick - prev_tick);
        min_ms = std::min(max_ms, current_tick - prev_tick);
        prev_tick = current_tick;

        // update the reported information
        if (current_tick - prev_report_tick > 1000)
        {
            reported_fps = frame_count * 1000. / double(current_tick - prev_report_tick);

            // 1/(ms/frame) = frame/ms
            // frame/ms * 1000ms/s = frame/s
            reported_min_fps = 1000. / min_ms;
            reported_max_fps = 1000. / max_ms;
            prev_report_tick = current_tick;
            max_ms = INT_MIN;
            min_ms = INT_MAX;
            frame_count = 0;
        }
    }

    [[nodiscard]] constexpr auto get_fps() const -> double { return reported_fps; }
    [[nodiscard]] constexpr auto get_min_fps() const -> double { return reported_min_fps; }
    [[nodiscard]] constexpr auto get_max_fps() const -> double { return reported_max_fps; }
};

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

/*
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
}*/

inline static constexpr glm::dvec3 CUBE_CORNER_1 = {-6, -6, -6};
inline static constexpr glm::dvec3 CUBE_CORNER_2 = {6, 6, 6};
inline static constexpr int64_t MARGIN_LEFT = 0;
inline static constexpr int64_t MARGIN_BOTTOM = 20;

class molorb_display
{
    using atoms_vertex_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{.enable_normal = true, .enable_color = true}>;
    using bound_box_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{.enable_normal = true, .enable_color = true, .primitive = GL_LINES}>;

    hartree_fock_result result;
    std::vector<hf_isosurface> iso;
    atoms_vertex_buffer atoms_vb;
    bound_box_buffer bounding_box_vb;
    size_t mo_index = 0;
    double rot_x = 0;
    double rot_y = 0;
    double zoom = 1;
    bool display_mo = true;
    bool display_atoms = true;
    bool display_box = true;

    void next_mo()
    {
        mo_index++;
        mo_index = (mo_index + result.orbitals.size()) % result.orbitals.size();
    }

    void prev_mo()
    {
        mo_index--;
        mo_index = (mo_index + result.orbitals.size()) % result.orbitals.size();
    }

    void generate_surfaces()
    {
        auto start = std::chrono::high_resolution_clock::now();
        iso.resize(result.orbitals.size());
        std::cout << "Generating your MOs\n";
        for (size_t i = 0; i < result.orbitals.size(); i++)
        {
            std::cout << "MO: " << i << '\n';
            iso[i].isolevel_hf(result, i, CUBE_CORNER_1, CUBE_CORNER_2, 20, 0.05);
        }
        auto stop = std::chrono::high_resolution_clock::now();

        std::cout << "Time taken by function: " << (duration_cast<std::chrono::microseconds>(stop - start)).count() << " microseconds" << std::endl;
    }

    void draw_sphere(float r, glm::vec3 center, uint32_t color)
    {
        static constexpr float di = 0.02;
        static constexpr float dj = 0.04;
        static constexpr float db = di * 2 * M_PI;
        static constexpr float da = dj * M_PI;

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
            draw_sphere(0.1, i.second, ELEMENT_COLORS[i.first - 1]);
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

    void draw_static(gl::render_manager& render)
    {
        if (!render.get_shader_manager().has_shader("atoms_shader"))
        {
            render.get_shader_manager()
                .add_shader("atoms_shader", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl")
                .compile_shader("atoms_shader");
            render.get_shader_manager()
                .push_shader("atoms_shader")
                .set_vec_float("diffuse_pos", {10, 10, 10})
                .set_vec_float("diffuse_color", {1, 1, 1})
                .set_vec_float("ambient_color", {1, 1, 1});
        }

        if (display_atoms)
        {
            render.get_shader_manager().push_shader("atoms_shader").set_float("ambient_strength", 0.2);
            render.render_vertex_buffer(atoms_vb);
            render.get_shader_manager().pop_shader();
        }

        if (display_box)
        {
            render.get_shader_manager().push_shader("atoms_shader").set_float("ambient_strength", 1);
            render.render_vertex_buffer(bounding_box_vb);
            render.get_shader_manager().pop_shader();
        }
    }

public:
    molorb_display(hartree_fock_result result) : result(std::move(result)) {}

    void init()
    {
        generate_atoms();
        generate_surfaces();
    }

    void key_pressed(int key, int scancode, int mods)
    {
        switch (key)
        {
        case GLFW_KEY_LEFT:
            prev_mo();
            break;
        case GLFW_KEY_RIGHT:
            next_mo();
            break;
        case GLFW_KEY_R:
            generate_surfaces();
            break;
        case GLFW_KEY_Z:
            display_mo = !display_mo;
            break;
        case GLFW_KEY_X:
            display_atoms = !display_atoms;
            break;
        case GLFW_KEY_C:
            display_box = !display_box;
            break;
        default:
            break;
        }
    }

    void mouse_dragged(double delta_x, double delta_y)
    {
        rot_x += delta_y * 0.005f;
        rot_y += delta_x * -0.005f;
    }

    void mouse_button_pressed(int button, int mods)
    {
        switch (button)
        {
        case GLFW_MOUSE_BUTTON_5:
            next_mo();
            break;
        case GLFW_MOUSE_BUTTON_4:
            prev_mo();
            break;
        default:
            break;
        }
    }

    void scrollwheel(double delta_x, double delta_y) { zoom = std::clamp(zoom + delta_y, -4., 4.); }

    void render(gl::render_manager& render)
    {
        render.view().push().mult(glm::lookAt(glm::vec3(10 * std::exp(zoom * 0.693), 0, 0), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0)));
        render.model().push().rotate(rot_x, 0, 0, 1).rotate(rot_y, 0, 1, 0);

        gl::gl_check();

        if (display_mo)
        {
            iso[mo_index].draw(render);
        }
        draw_static(render);

        render.model().pop();
        render.view().pop();

        static constexpr auto EV_PER_HT = 27.211407953;

        int electron_in_orbital = int(result.electron_count) - int(mo_index * 2);
        std::string electron_desc;

        if (mo_index == result.homo_index)
        {
            electron_desc = fmt::format("Electrons: {} (HOMO)", electron_in_orbital);
        }
        else if (mo_index == result.homo_index + 1)
        {
            electron_desc = "Electrons: 0 (LUMO)";
        }
        else if (electron_in_orbital > 0)
        {
            electron_desc = fmt::format("Electrons: {} (Occupied, HOMO-{})", std::min(electron_in_orbital, 2), result.homo_index - mo_index);
        }
        else
        {
            electron_desc = fmt::format("Electrons: 0 (Unoccupied, LUMO+{})", mo_index - result.homo_index - 1);
        }

        gl::draw_text(fmt::format(
            "Info:\n"
            "  Iterations: {}\n"
            "  Electrons: {}\n"
            "MO Info:\n"
            "  MO #{}\n"
            "  E: {:.4f} eV\n"
            "  {}", result.iterations,result.electron_count, mo_index + 1, result.mo_energies[mo_index] * EV_PER_HT, electron_desc), {3, 0});
    }
};

static auto format_status(const fps_calculator& fps) -> std::string
{
    return fmt::format(PROJECT_NAME " v" VERSION " - c++ = " MESON_CXX_COMPILER " - {:.2f} MiB - {}/{}/{} FPS", memory_used() / double(1024),
                       (size_t)fps.get_fps(), (size_t)fps.get_max_fps(), (size_t)fps.get_min_fps());
}

class mol_orbital : public gl::window
{
    gl::render_manager man;
    molorb_display molorb_ui;
    bool is_left_pressed = false;

    double prev_mouse_x = NAN;
    double prev_mouse_y{};
    fps_calculator fps;

protected:
    void init() override
    {
        // glDebugMessageCallback(debugMessageCallback, nullptr);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
        gl::gl_check();

        fov = 5;
        z_near = 0.1;
        z_far = 100;
        bg = {0.1, 0.1, 0.1, 1};
        man.get_shader_manager().add_shader("main", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl").compile_shader("main");

        molorb_ui.init();
    }

    void key_pressed(int key, int scancode, int mods) override
    {
        switch (key)
        {
        case GLFW_KEY_ESCAPE:
            close();
            break;
        default:
            molorb_ui.key_pressed(key, scancode, mods);
        }
    }

    void key_released(int key, int scancode, int mods) override {}

    void mouse_moved(double xpos, double ypos) override
    {
        double dx = xpos - prev_mouse_x;
        double dy = ypos - prev_mouse_y;

        if (dx != NAN)
        {
            if (is_left_pressed)
            {
                molorb_ui.mouse_dragged(dx, dy);
            }
        }

        prev_mouse_x = xpos;
        prev_mouse_y = ypos;
    }

    void mouse_button_pressed(int button, int mods) override
    {
        switch (button)
        {
        case GLFW_MOUSE_BUTTON_LEFT:
            is_left_pressed = true;
            break;
        default:
            molorb_ui.mouse_button_pressed(button, mods);
        }
    }

    void mouse_button_released(int button, int mods) override
    {
        switch (button)
        {
        case GLFW_MOUSE_BUTTON_LEFT:
            is_left_pressed = false;
            break;
        default:
            break;
        }
    }

    void mouse_scrolled(double x_off, double y_off) override { molorb_ui.scrollwheel(x_off, y_off); }

    void render(const gl::mat4& proj, glm::uvec2 screen_size) override
    {
        glViewport(0, 0, screen_size.x, screen_size.y);
        man.get_shader_manager().use_shader("main");
        man.projection().push().mult(proj);
        man.set_screen_size(screen_size);
        molorb_ui.render(man);
        man.projection().pop();

        gl::draw_text(format_status(fps), {0, screen_size.y - MARGIN_BOTTOM});
        fps.report_frame();
    }

public:
    mol_orbital(int a, int b, const char* c, hartree_fock_result result) : gl::window(a, b, c), molorb_ui(std::move(result)) {}
};

auto main(int argc, char* argv[]) -> int
{
    try
    {
        std::string filename;

        if (argc != 2)
        {
            std::cerr << "usage: " << argv[0] << " [filename]\n";
            exit(-1);
        }

        filename = argv[1];

        basis_manager::get_instance().register_basis("sto-2g");
        basis_manager::get_instance().register_basis("sto-3g");
        basis_manager::get_instance().register_basis("sto-4g");
        basis_manager::get_instance().register_basis("sto-5g");
        basis_manager::get_instance().register_basis("sto-6g");

        auto mol = molecule::read_from_file(filename);

        std::cout << "solving\n";
        auto res = solve(mol);
        mol_orbital window(1000, 1000, "Molecular Orbitals", std::move(res));

        window.run();
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

