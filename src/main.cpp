#include <GL/glew.h>
//
#include "build_config.h"
#include "fps.h"
#include "gl/2dui.h"
#include "gl/2dui_builder.h"
#include "gl/event.h"
#include "gl/font.h"
#include "gl/orbital_render.h"
#include "gl/vertex_buffer.h"
#include "gl/window.h"
#include "hf/hf.h"
#include "keybind.h"
#include "nfd.h"
#include "platform_dep.h"
#include "solve_in_thread.h"
#include "ui/config_ui.h"
#include "ui/electron_cloud_manager.h"
#include "ui/molorb_display.h"
#include <GLFW/glfw3.h>
#include <chrono>
#include <filesystem>
#include <fmt/core.h>
#include <glm/glm.hpp>
#include <memory>
#include <nlohmann/json.hpp>

/*
inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_SOURCE[] = {{GL_DEBUG_SOURCE_API, "api"},
                                                                            {GL_DEBUG_SOURCE_WINDOW_SYSTEM, "window system"},
                                                                            {GL_DEBUG_SOURCE_SHADER_COMPILER, "shader compiler"},
                                                                            {GL_DEBUG_SOURCE_THIRD_PARTY, "third party"},
                                                                            {GL_DEBUG_SOURCE_APPLICATION, "application"},
                                                                            {GL_DEBUG_SOURCE_OTHER, "other"}};

inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_TYPE[] = {{GL_DEBUG_TYPE_ERROR, "error"},
                                                                          {GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR, "deprecated behavior"},
                                                                          {GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR, "undefined behavior"},
                                                                          {GL_DEBUG_TYPE_PORTABILITY, "portability"},
                                                                          {GL_DEBUG_TYPE_PERFORMANCE, "performance"},
                                                                          {GL_DEBUG_TYPE_MARKER, "marker"},
                                                                          {GL_DEBUG_TYPE_PUSH_GROUP, "push group"},
                                                                          {GL_DEBUG_TYPE_POP_GROUP, "pop group"},
                                                                          {GL_DEBUG_TYPE_OTHER, "other"}};

inline static constexpr std::pair<GLenum, const char*> GL_DEBUG_LEVEL[] = {{GL_DEBUG_SEVERITY_HIGH, "high"},
                                                                           {GL_DEBUG_SEVERITY_MEDIUM, "medium"},
                                                                           {GL_DEBUG_SEVERITY_LOW, "low"},
                                                                           {GL_DEBUG_SEVERITY_NOTIFICATION, "notification"}};

auto string_lookup(GLenum value, const auto& stringList) -> const char*
{
    for (const auto& pair : stringList)
    {
        if (pair.first == value)
        {
            return pair.second;
        }
    }
    return "unknown";
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
}*/

inline static constexpr auto MiB_SIZE = 1024 * 1024;

namespace
{
    auto format_status(const fps_calculator& fps) -> std::string
    {
        return fmt::format(PROJECT_NAME " v" VERSION " - " MESON_CXX_COMPILER "/" MESON_C_COMPILER " - {:.2f} MiB - {}/{}/{} FPS",
                           memory_used() / double(MiB_SIZE), (size_t)fps.get_fps(), (size_t)fps.get_max_fps(), (size_t)fps.get_min_fps());
    }
} // namespace

class mol_orbital : public gl::window
{
    gl::render_manager renderer;
    molorb_display molorb_ui;
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
        renderer.get_shader_manager().add_shader("main", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl").compile_shader("main");

        molorb_ui.init();
    }

    void key_pressed(int key, int scancode, int mods) override
    {
        if (molorb_ui.handle_key_press(key, scancode, mods))
        {
            return;
        }

        switch (key)
        {
        case GLFW_KEY_ESCAPE:
            close();
            break;
        }
    }

    void key_released(int key, int scancode, int mods) override { molorb_ui.handle_key_release(key, scancode, mods); }

    void mouse_moved(double mouse_x, double mouse_y) override { molorb_ui.handle_mouse(mouse_x, mouse_y); }

    void mouse_button_pressed(int button, int mods) override { molorb_ui.handle_mouse_press(button, mods); }

    void mouse_button_released(int button, int mods) override { molorb_ui.handle_mouse_release(button, mods); }

    void mouse_scrolled(double x_off, double y_off) override { molorb_ui.scrollwheel(x_off, y_off); }

    void render(const gl::mat4& proj, glm::uvec2 screen_size) override
    {
        glViewport(0, 0, screen_size.x, screen_size.y);
        renderer.get_shader_manager().use_shader("main");
        renderer.projection().push().mult(proj);
        renderer.set_screen_size(screen_size);
        molorb_ui.render(renderer);
        renderer.projection().pop();

        gl::render_text(format_status(fps), {3, screen_size.y - MARGIN_BOTTOM});
        fps.report_frame();
    }

public:
    mol_orbital(int width, int height, const char* name,
                std::variant<hartree_fock_result, std::future<hartree_fock_result>, std::monostate> init_state)
        : gl::window(width, height, name), molorb_ui(std::move(init_state))
    {
    }
};

auto main(int argc, const char* argv[]) -> int
{
    std::span<const char*> args(argv, argc);
    try
    {
        std::filesystem::create_directories(get_home() + "/.mo_config/");
        keybind_config::get_instance().load(get_home() + "/.mo_config/keybind.json");

        basis_manager::get_instance().register_basis("sto-2g");
        basis_manager::get_instance().register_basis("sto-3g");
        basis_manager::get_instance().register_basis("sto-4g");
        basis_manager::get_instance().register_basis("sto-5g");
        basis_manager::get_instance().register_basis("sto-6g");

        if (argc == 1)
        {
            mol_orbital window(1000, 1000, "Molecular Orbitals", std::monostate{});
            window.run();
        }
        else if (argc == 2)
        {
            mol_orbital window(1000, 1000, "Molecular Orbitals", run_hf(args[1]));
            window.run();
        }
        else
        {
            std::cerr << "usage: " << args[0] << " [filename]\n";
            exit(-1);
        }
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
}

