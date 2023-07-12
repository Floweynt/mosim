#include "ui/molorb_display.h"
#include "keybind.h"
#include "nfd.h"
#include "ui_config.h"

inline static constexpr auto EV_PER_HT = 27.211407953;

void molorb_display::bake_sphere(float radius, glm::vec3 center, uint32_t color)
{
    static constexpr float DELTA_B = SPHERE_SUBDIV_A * 2 * M_PI;
    static constexpr float DELTA_A = SPHERE_SUBDIV_B * M_PI;

    for (float i = 0; i < 1.0; i += SPHERE_SUBDIV_A)
    {
        for (float j = 0; j < 1.0; j += SPHERE_SUBDIV_B)
        {
            float b = i * 2 * M_PI;
            float a = (j - 0.5) * M_PI;

            glm::vec3 normal = {std::cos(a + DELTA_A / 2) * std::cos(b + DELTA_B / 2), std::cos(a + DELTA_A / 2) * std::sin(b + DELTA_B / 2),
                                std::sin(a + DELTA_A / 2)};
            glm::vec3 pt0(radius * std::cos(a) * std::cos(b), radius * std::cos(a) * std::sin(b), radius * std::sin(a));
            glm::vec3 pt1(radius * std::cos(a) * std::cos(b + DELTA_B), radius * std::cos(a) * std::sin(b + DELTA_B), radius * std::sin(a));
            glm::vec3 pt2(radius * std::cos(a + DELTA_A) * std::cos(b + DELTA_B), radius * std::cos(a + DELTA_A) * std::sin(b + DELTA_B),
                          radius * std::sin(a + DELTA_A));
            glm::vec3 pt3(radius * std::cos(a + DELTA_A) * std::cos(b), radius * std::cos(a + DELTA_A) * std::sin(b), radius * std::sin(a + DELTA_A));

            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt0 + center).color(color).normal(normal).end());
            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt1 + center).color(color).normal(normal).end());
            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt2 + center).color(color).normal(normal).end());

            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt2 + center).color(color).normal(normal).end());
            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt3 + center).color(color).normal(normal).end());
            atoms_vb.vert(atoms_vertex_buffer::builder().vert(pt0 + center).color(color).normal(normal).end());
        }
    }
}

#define select(n, v) (edge & (1 << (n))) ? (CUBE_CORNER_1.v) : (CUBE_CORNER_2.v)

void molorb_display::generate_static()
{
    // reset the vb
    atoms_vb = atoms_vertex_buffer();
    for (const auto& atom : result.atoms)
    {
        bake_sphere(ATOM_SPHERE_RADIUS, atom.second, ELEMENT_COLORS[atom.first - 1]);
    }
    atoms_vb.bake();

    static constexpr std::array<uint8_t, 12> EDGES = {0b000100, 0b100110, 0b110010, 0b010000, 0b001101, 0b101111,
                                                      0b111011, 0b011001, 0b000001, 0b100101, 0b110111, 0b010011};

    for (auto edge : EDGES)
    {
        bounding_box_vb.vert(bound_box_buffer::builder().vert(select(0, x), select(1, y), select(2, z)).color(BOX_COLOR).end());
        bounding_box_vb.vert(bound_box_buffer::builder().vert(select(3, x), select(4, y), select(5, z)).color(BOX_COLOR).end());
    }
    bounding_box_vb.bake();
}

void molorb_display::render_static(gl::render_manager& render)
{
    if (!render.get_shader_manager().has_shader("static_shader"))
    {
        render.get_shader_manager()
            .add_shader("static_shader", "assets/shaders/vertex.glsl", "assets/shaders/fragment.glsl")
            .compile_shader("static_shader");
        render.get_shader_manager()
            .push_shader("static_shader")
            .set_vec_float("diffuse_pos", LIGHT_POS)
            .set_vec_float("diffuse_color", {TEXT_COLOR})
            .set_vec_float("ambient_color", {TEXT_COLOR});
    }

    if (display_atoms)
    {
        render.get_shader_manager().push_shader("static_shader").set_float("ambient_strength", 0.2);
        render.render_vertex_buffer(atoms_vb);
        render.get_shader_manager().pop_shader();
    }

    if (display_box)
    {
        render.get_shader_manager().push_shader("static_shader").set_float("ambient_strength", 1);
        render.render_vertex_buffer(bounding_box_vb);
        render.get_shader_manager().pop_shader();
    }
}

// EVENT HANDLERS
auto molorb_display::key_pressed(int key, int /*scancode*/, int mods) -> bool
{
    auto& kb = keybind_config::get_instance();
    bool flag = false;
    flag |= on_keybind(kb.goto_lumo, key, mods, [this]() { electron_cloud.set_mo(result.homo_index + 1);}, false);
    flag |= on_keybind(kb.goto_homo, key, mods, [this]() { electron_cloud.set_mo(result.homo_index);}, false);
    flag |= on_keybind(kb.goto_next_mo, key, mods, [this]() { electron_cloud.next_mo();}, false);
    flag |= on_keybind(kb.goto_prev_mo, key, mods, [this]() { electron_cloud.prev_mo();}, false);
    flag |= on_keybind(kb.goto_next_mo_alt, key, mods, [this]() { electron_cloud.next_mo();}, false);
    flag |= on_keybind(kb.goto_prev_mo_alt, key, mods, [this]() { electron_cloud.prev_mo();}, false);
    flag |= on_keybind(kb.open_config, key, mods, [this]() { config->set_active(true); }, false);
    flag |= on_keybind(kb.rerender, key, mods, [this]() { electron_cloud.regenerate();}, false);
    flag |= on_keybind(kb.toggle_render_mo, key, mods, [this]() { display_mo = !display_mo; }, false);
    flag |= on_keybind(kb.toggle_render_atoms, key, mods, [this]() { display_atoms = !display_atoms; }, false);
    flag |= on_keybind(kb.toggle_render_box, key, mods, [this]() { display_box = !display_box; }, false);
    flag |= on_keybind(kb.save_solution, key, mods, [this]() { save(); }, false);
    return flag;
}

auto molorb_display::mouse_dragged(double delta_x, double delta_y) -> bool
{
    rot_x += delta_y * MOUSE_MOVE_FACTOR;
    rot_y += delta_x * -MOUSE_MOVE_FACTOR;
    return true;
}

auto molorb_display::mouse_button_pressed(int button, int mods) -> bool
{
    auto& kb = keybind_config::get_instance();
    bool flag = false;
    flag |= on_keybind(kb.goto_lumo, button, mods, [this]() { electron_cloud.set_mo(result.homo_index + 1);}, true);
    flag |= on_keybind(kb.goto_homo, button, mods, [this]() { electron_cloud.set_mo(result.homo_index);}, true);
    flag |= on_keybind(kb.goto_next_mo, button, mods, [this]() { electron_cloud.next_mo();}, true);
    flag |= on_keybind(kb.goto_prev_mo, button, mods, [this]() { electron_cloud.prev_mo();}, true);
    flag |= on_keybind(kb.goto_next_mo_alt, button, mods, [this]() { electron_cloud.next_mo();}, true);
    flag |= on_keybind(kb.goto_prev_mo_alt, button, mods, [this]() { electron_cloud.prev_mo();}, true);
    flag |= on_keybind(kb.open_config, button, mods, [this]() { config->set_active(true); }, true);
    flag |= on_keybind(kb.rerender, button, mods, [this]() { electron_cloud.regenerate();}, true);
    flag |= on_keybind(kb.toggle_render_mo, button, mods, [this]() { display_mo = !display_mo; }, true);
    flag |= on_keybind(kb.toggle_render_atoms, button, mods, [this]() { display_atoms = !display_atoms; }, true);
    flag |= on_keybind(kb.toggle_render_box, button, mods, [this]() { display_box = !display_box; }, true);
    flag |= on_keybind(kb.save_solution, button, mods, [this]() { save(); }, true);
    return flag;
}

molorb_display::molorb_display(hartree_fock_result result) : result(std::move(result)), config(std::make_shared<config_ui>("a"))
{
    mo_display_buttons = gl::build_glue(
        "mo_display_buttons_glue", gl::ui_glue::GLUE_BOTTOM_LEFT,
        gl::padding_builder("mo_display_buttons_padding", gl::vertical_stack_builder("mo_display_buttons")
                                                              .margin(5)
                                                              .add_child(gl::button_builder("goto_homo", "Goto HOMO")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                     electron_cloud.set_mo(this->result.homo_index);
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                                              .add_child(gl::button_builder("goto_lumo", "Goto LUMO")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                     electron_cloud.set_mo(this->result.homo_index + 1);
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                                              .add_child(gl::button_builder("goto_next_mo", "Next MO")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                     electron_cloud.next_mo();
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                                              .add_child(gl::button_builder("goto_prev_mo", "Prev MO")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                     electron_cloud.prev_mo();
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                                              .add_child(gl::button_builder("rerender", "Rerender")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                     electron_cloud.regenerate();
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                                              .build())
            .pad(3)
            .bottom(MARGIN_BOTTOM + 3)
            .build());

    config_button = gl::build_glue("config_glue", gl::ui_glue::GLUE_TOP_RIGHT,
                                   gl::padding_builder("config_padding", gl::button_builder("config", "Open Config")
                                                                             .size(100, 20)
                                                                             .bg_color(PASSIVE_COLOR)
                                                                             .hover_color(ACTIVE_COLOR)
                                                                             .text_color(TEXT_COLOR)
                                                                             .on_press([this](int key, int /*mods*/) {
                                                                                 if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                                 {
                                                                                 config->set_active(true);
                                                                                     return true;
                                                                                 }
                                                                                 return false;
                                                                             })
                                                                             .build())
                                       .pad(3)
                                       .build());

    drag_handler = std::make_shared<gl::mouse_drag_handler<GLFW_MOUSE_BUTTON_LEFT>>();

    add_child(mo_display_buttons);
    add_child(config);
    add_child(config_button);
    add_child(drag_handler);
    set_key_press_handler([this](int key, int scancode, int mods) { return key_pressed(key, scancode, mods); });
    set_mouse_press_handler([this](int key, int mods) { return mouse_button_pressed(key, mods); });
    drag_handler->set_handler([this](float delta_x, float delta_y) { return mouse_dragged(delta_x, delta_y); });
}

void molorb_display::init()
{
    generate_static();
    electron_cloud.set_result(result);
}

void molorb_display::save()
{
    char* out_path = nullptr;
    nfdfilteritem_t filter[1] = {{"Molecular Orbital Output File", "json"}};
    if (NFD_SaveDialog(&out_path, filter, 1, ".", "mo_solution.json") != NFD_OKAY)
    {
        return;
    }

    std::string out(out_path);
    NFD_FreePathN(out_path);
    std::ofstream ofs(out);
    if (!ofs)
    {
        throw std::runtime_error("failed to save to file " + out);
    }

    nlohmann::json json;
    write_result(result, json);
    ofs << json;
}

void molorb_display::scrollwheel(double delta_x, double delta_y) { zoom = std::clamp(zoom + delta_y, -4., 4.); }

[[nodiscard]] auto molorb_display::get_electron_desc() const -> std::string
{
    int electron_in_orbital = int(result.electron_count) - int(electron_cloud.get_index() * 2);
    std::string electron_desc;

    if (electron_cloud.get_index() == result.homo_index)
    {
        electron_desc = fmt::format("Electrons: {} (HOMO)", electron_in_orbital);
    }
    else if (electron_cloud.get_index() == result.homo_index + 1)
    {
        electron_desc = "Electrons: 0 (LUMO)";
    }
    else if (electron_in_orbital > 0)
    {
        electron_desc =
            fmt::format("Electrons: {} (Occupied, HOMO-{})", std::min(electron_in_orbital, 2), result.homo_index - electron_cloud.get_index());
    }
    else
    {
        electron_desc = fmt::format("Electrons: 0 (Unoccupied, LUMO+{})", electron_cloud.get_index() - result.homo_index - 1);
    }

    return electron_desc;
}

void molorb_display::render(gl::render_manager& renderer)
{
    renderer.view().push().mult(glm::lookAt(glm::vec3(10 * std::exp(zoom * MOUSE_ZOOM_FACTOR), 0, 0), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0)));
    renderer.model().push().rotate(rot_x, 0, 0, 1).rotate(rot_y, 0, 1, 0);

    gl::gl_check();

    if (display_mo)
    {
        electron_cloud.render(renderer);
    }

    render_static(renderer);

    renderer.model().pop();
    renderer.view().pop();

    gl::render_text(fmt::format("Info:\n"
                                "  Iterations: {}\n"
                                "  Electrons: {}\n"
                                "MO Info:\n"
                                "  MO #{}\n"
                                "  E: {:.4f} eV\n"
                                "  {}",
                                result.iterations, result.electron_count, electron_cloud.get_index() + 1,
                                result.mo_energies[electron_cloud.get_index()] * EV_PER_HT, get_electron_desc()),
                    {3, 3});

    mo_display_buttons->render(renderer);
    config_button->render(renderer);
    config->render(renderer);
}

