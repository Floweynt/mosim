#include "ui/config_ui.h"

inline static constexpr std::pair<int, const char*> KEY_NAMES[] = {
    {GLFW_KEY_SPACE, "Space"},
    {GLFW_KEY_GRAVE_ACCENT, "Grave"},
    {GLFW_KEY_ESCAPE, "Esc"},
    {GLFW_KEY_ENTER, "Enter"},
    {GLFW_KEY_TAB, "Tab"},
    {GLFW_KEY_BACKSPACE, "Bsp"},
    {GLFW_KEY_INSERT, "Ins"},
    {GLFW_KEY_DELETE, "Del"},
    {GLFW_KEY_RIGHT, "Right"},
    {GLFW_KEY_LEFT, "Left"},
    {GLFW_KEY_DOWN, "Down"},
    {GLFW_KEY_UP, "Up"},
    {GLFW_KEY_PAGE_UP, "PgUp"},
    {GLFW_KEY_PAGE_DOWN, "PgDn"},
    {GLFW_KEY_HOME, "Home"},
    {GLFW_KEY_END, "End"},
    {GLFW_KEY_CAPS_LOCK, "caps_lock"},
    {GLFW_KEY_SCROLL_LOCK, "scroll_lock"},
    {GLFW_KEY_NUM_LOCK, "num_lock"},
    {GLFW_KEY_PRINT_SCREEN, "PrnSc"},
    {GLFW_KEY_PAUSE, "Pause"},
};

inline static constexpr std::pair<int, const char*> MOUSE_BUTTON_NAMES[] = {
    {GLFW_MOUSE_BUTTON_LEFT, "Left-c"},
    {GLFW_MOUSE_BUTTON_RIGHT, "Right-c"},
    {GLFW_MOUSE_BUTTON_MIDDLE, "Mid-c"},
    {GLFW_MOUSE_BUTTON_4, "Btn4"},
    {GLFW_MOUSE_BUTTON_5, "Btn5"},
    {GLFW_MOUSE_BUTTON_6, "Btn6"},
    {GLFW_MOUSE_BUTTON_7, "Btn7"},
    {GLFW_MOUSE_BUTTON_8, "Btn8"},
};

inline static constexpr std::pair<int, const char*> MOD_NAMES[] = {
    {GLFW_MOD_SUPER, "Su-"},
    {GLFW_MOD_CONTROL, "C-"},
    {GLFW_MOD_ALT, "A-"},
    {GLFW_MOD_SHIFT, "S-"},
};

static auto kb_to_string(const keybind& kb) -> std::string
{
    if (kb.key == -1)
    {
        return "[none]";
    }

    std::string kb_id;
    for (const auto& entry : MOD_NAMES)
    {
        if (kb.mods & entry.first)
        {
            kb_id += entry.second;
        }
    }

    if (kb.is_mouse)
    {
        const char* name = "UNK";
        for (const auto& entry : MOUSE_BUTTON_NAMES)
        {
            if (entry.first == kb.key)
            {
                name = entry.second;
                break;
            }
        }

        return kb_id + name;
    }

    const auto* name_cstr = glfwGetKeyName(kb.key, kb.scan);
    std::string name = name_cstr != nullptr ? glfwGetKeyName(kb.key, kb.scan) : "";
    if (name.empty())
    {

        for (const auto& entry : KEY_NAMES)
        {
            if (entry.first == kb.key)
            {
                name = entry.second;
                break;
            }
        }
    }

    if (name.empty())
    {
        if (kb.key >= GLFW_KEY_F1 && kb.key <= GLFW_KEY_F25)
        {
            name = fmt::format("F{}", kb.key - GLFW_KEY_F1 + 1);
        }
    }

    if (name.empty())
    {
        return fmt::format("{}-{}-{}", kb.key, kb.scan, kb.mods);
    }

    return kb_id + name;
}

void keybind_ui::set_keybind_target(keybind* val)
{
    if (keybind_target != val)
    {
        keybind_target = val;
        mark_dirty();
    }
}

void keybind_ui::set_hovering_target(keybind* val)
{
    if (hovering_target != val)
    {
        hovering_target = val;
        mark_dirty();
    }
}

void keybind_ui::update(glm::uvec2 win_size)
{
    // the ui looks like this
    // +-----------------------------------+
    // |             +--------+ +-------+  |
    // |  label      | MainKb | | AltKb |  |
    // |             +--------+ +-------+  |
    // *-----------------------------------+

    draw_rect({0, 0}, KB_UI_SIZE, {0.3, 0.3, 0.3});

    auto middle_y = (get_size() / 2.f).y;
    draw_text({LABEL_OFFSET, middle_y}, label, TEXT_COLOR, 1.1, TEXT_LEFT, TEXT_CENTER);
    draw_rect(MAIN_KB_POS_1, MAIN_KB_POS_2, hovering_target == main ? ACTIVE_COLOR : PASSIVE_COLOR);
    draw_text((MAIN_KB_POS_1 + MAIN_KB_POS_2) / 2.f, keybind_target == main ? "?" : kb_to_string(*main), TEXT_COLOR, 1, TEXT_CENTER, TEXT_CENTER);
    if (alt != nullptr)
    {
        draw_rect(ALT_KB_POS_1, ALT_KB_POS_2, hovering_target == alt ? ACTIVE_COLOR : PASSIVE_COLOR);
        draw_text((ALT_KB_POS_1 + ALT_KB_POS_2) / 2.f, keybind_target == alt ? "?" : kb_to_string(*alt), TEXT_COLOR, 1, TEXT_CENTER, TEXT_CENTER);
    }
}

auto keybind_ui::handle_mouse(float x, float y) -> bool
{
    glm::vec2 mouse_pos(x, y);

    if (in_bb_corners(pos() + MAIN_KB_POS_1, pos() + MAIN_KB_POS_2, mouse_pos))
    {
        set_hovering_target(main);
        return true;
    }

    if (in_bb_corners(pos() + ALT_KB_POS_1, pos() + ALT_KB_POS_2, mouse_pos) && (alt != nullptr))
    {
        set_hovering_target(alt);
        return true;
    }

    set_hovering_target(nullptr);
    return false;
}

auto keybind_ui::handle_mouse_press(int key, int mods) -> bool
{
    if (keybind_target != nullptr)
    {
        *keybind_target = keybind{true, key, 0, mods & GLFW_MASK_RELEVENT};
        on_update();
        set_keybind_target(nullptr);
        return true;
    }

    if (key == GLFW_MOUSE_BUTTON_LEFT && (hovering_target != nullptr))
    {
        set_keybind_target(hovering_target);
        return true;
    }

    return false;
}

auto keybind_ui::handle_key_press(int key, int scancode, int mods) -> bool
{
    if (keybind_target != nullptr)
    {
        if (key >= GLFW_KEY_LEFT_SHIFT || (key >= GLFW_KEY_CAPS_LOCK && key <= GLFW_KEY_NUM_LOCK))
        {
            return true;
        }

        *keybind_target = keybind{false, key, scancode, mods};
        on_update();
        set_keybind_target(nullptr);
        return true;
    }

    return false;
}

keybind_ui::keybind_ui(std::string name, std::string label, std::function<void()> on_update, keybind* main, keybind* alt)
    : ui_widget("keybind", std::move(name)), label(std::move(label)), main(main), alt(alt), on_update(std::move(on_update))
{
    set_size(KB_UI_SIZE);
}

void config_ui::update(glm::uvec2 win_size)
{
    if (!is_active)
    {
        return;
    }
    keybind_stack->move(pos());
    update_children(win_size);
    set_size(keybind_stack->get_size());
    move(glm::vec2(win_size) / 2.f - get_size() / 2.f);
    draw_rect({0, 0}, get_size(), CONFIG_BG_COLOR);
    draw_text({CONFIG_UI_SIZE.x / 2, 3}, "Config Menu", TEXT_COLOR, 1.2, TEXT_CENTER, TEXT_TOP);
}

void config_ui::render_hook(gl::render_manager& renderer)
{
    if (!is_active)
    {
        return;
    }

    for (const auto& child : get_children())
    {
        child->render(renderer);
    }
}

void config_ui::on_update() { keybind_config::get_instance().save(get_home() + "/.mo_config/keybind.json"); }

config_ui::config_ui(std::string name) : ui_widget("config_menu", std::move(name))
{
    auto& kb = keybind_config::get_instance();
    keybind_stack = gl::padding_builder("config_menu_padding",
        gl::vertical_stack_builder("config_menu_stack")
            .margin(1)
            .add_child(gl::padding_builder("close_button_padding", gl::button_builder("close_button", "x")
                                                                       .bg_color(PASSIVE_COLOR)
                                                                       .hover_color(ACTIVE_COLOR)
                                                                       .text_color(TEXT_COLOR)
                                                                       .size(20, 20)
                                                                       .on_press([this](int key, int) {
                                                                           if (key == GLFW_MOUSE_BUTTON_LEFT)
                                                                           {
                                                                               set_active(false);
                                                                               return true;
                                                                           }
                                                                           return false;
                                                                       })
                                                                       .build())
                           .bottom(3)
                           .build())
            .add_child(
                gl::label_builder("keybind_label", "Keybinds").size(CONFIG_UI_SIZE.x, 30).bg_color(CONFIG_BG_COLOR).text_color(TEXT_COLOR).build())
            .add_child(std::make_shared<keybind_ui>("keybind_homo", "Goto HOMO", on_update, &kb.goto_homo))
            .add_child(std::make_shared<keybind_ui>("keybind_lumo", "Goto HOMO", on_update, &kb.goto_lumo))
            .add_child(std::make_shared<keybind_ui>("keybind_rerender", "Rerender", on_update, &kb.rerender))
            .add_child(std::make_shared<keybind_ui>("keybind_next_mo", "Goto next MO", on_update, &kb.goto_next_mo, &kb.goto_next_mo_alt))
            .add_child(std::make_shared<keybind_ui>("keybind_prev_mo", "Goto previous MO", on_update, &kb.goto_prev_mo, &kb.goto_prev_mo_alt))
            .add_child(std::make_shared<keybind_ui>("keybind_open_config", "Open Config", on_update, &kb.open_config))
            .add_child(std::make_shared<keybind_ui>("keybind_toggle_render_mo", "Toggle render MO", on_update, &kb.toggle_render_mo))
            .add_child(std::make_shared<keybind_ui>("keybind_toggle_render_atoms", "Toggle render atoms", on_update, &kb.toggle_render_atoms))
            .add_child(std::make_shared<keybind_ui>("keybind_toggle_render_box", "Toggle render box", on_update, &kb.toggle_render_box))
            .add_child(std::make_shared<keybind_ui>("keybind_save", "Save MO coefficients", on_update, &kb.save_solution))
            .build())
                        .pad(3)
                        .build();
    add_child(keybind_stack);
}

auto config_ui::handle_mouse(float x, float y) -> bool
{
    auto ptr = std::dynamic_pointer_cast<gl::event_handler>(keybind_stack);
    ptr->handle_mouse(x, y);
    return is_active;
}

auto config_ui::handle_mouse_press(int key, int mods) -> bool
{
    if (!is_active)
    {
        return false;
    }
    auto ptr = std::dynamic_pointer_cast<gl::event_handler>(keybind_stack);
    ptr->handle_mouse_press(key, mods);

    return true;
}

auto config_ui::handle_mouse_release(int key, int mods) -> bool
{
    if (!is_active)
    {
        return false;
    }
    auto ptr = std::dynamic_pointer_cast<gl::event_handler>(keybind_stack);
    ptr->handle_mouse_release(key, mods);

    return true;
}

auto config_ui::handle_key_press(int key, int scancode, int mods) -> bool
{
    if (!is_active)
    {
        return false;
    }

    auto ptr = std::dynamic_pointer_cast<gl::event_handler>(keybind_stack);
    if (ptr->handle_key_press(key, scancode, mods))
    {
        return true;
    }

    if (key == GLFW_KEY_ESCAPE)
    {
        set_active(false);
    }

    return true;
}

auto config_ui::handle_key_release(int key, int scancode, int mods) -> bool
{
    if (!is_active)
    {
        return false;
    }
    auto ptr = std::dynamic_pointer_cast<gl::event_handler>(keybind_stack);
    ptr->handle_key_release(key, scancode, mods);
    return true;
}

void config_ui::set_active(bool flag)
{
    if (is_active != flag)
    {
        mark_dirty();
        is_active = flag;
    }
}
