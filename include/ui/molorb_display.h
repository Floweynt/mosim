#pragma once


#include "gl/2dui.h"
#include "gl/2dui_builder.h"
#include "gl/event.h"
#include "gl/vertex_buffer.h"
#include "hf/hf.h"
#include "ui/config_ui.h"
#include "ui_config.h"
#include "ui/electron_cloud_manager.h"
#include <GLFW/glfw3.h>
#include <future>
#include <variant>

class molorb_display : public gl::compound_event_handler
{
    using atoms_vertex_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{.enable_normal = true, .enable_color = true}>;
    using bound_box_buffer = gl::vertex_buffer<gl::vertex_buffer_cfg{.enable_normal = true, .enable_color = true, .primitive = GL_LINES}>;

    std::variant<hartree_fock_result, std::future<hartree_fock_result>, std::monostate> result;
    electron_cloud_manager electron_cloud;
    atoms_vertex_buffer atoms_vb;
    bound_box_buffer bounding_box_vb;
    double rot_x = 0;
    double rot_y = 0;
    double zoom = 1;
    bool display_mo = true;
    bool display_atoms = true;
    bool display_box = true;
    std::shared_ptr<config_ui> config;

    // UI components
    std::shared_ptr<gl::mouse_drag_handler<GLFW_MOUSE_BUTTON_LEFT>> drag_handler;
    gl::widget_ref mo_display_buttons;
    gl::widget_ref config_button;
    gl::widget_ref file_ops_button;

    void bake_sphere(float radius, glm::vec3 center, uint32_t color);
    void generate_static();
    void render_static(gl::render_manager& render);
    // EVENT HANDLERS
    auto key_pressed(int key, int /*scancode*/, int mods) -> bool;
    auto mouse_dragged(double delta_x, double delta_y) -> bool;
    auto mouse_button_pressed(int button, int /*mods*/) -> bool;
    void save();
    void load();
    void reload();
public:
    molorb_display(std::variant<hartree_fock_result, std::future<hartree_fock_result>, std::monostate> init_state);
    void init();
    void scrollwheel(double delta_x, double delta_y);
    [[nodiscard]] auto get_electron_desc() const -> std::string;
    void render(gl::render_manager& renderer);
};
