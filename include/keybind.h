#pragma once

#include "singleton.h"
#include <GLFW/glfw3.h>
#include <fstream>
#include <functional>
#include <ios>
#include <string>

struct keybind
{
    bool is_mouse{};
    int key = -1;
    int scan{};
    int mods{};
};

class keybind_config : public singleton<keybind_config>
{
public:
    keybind goto_lumo{.key = GLFW_KEY_L};
    keybind goto_homo{.key = GLFW_KEY_H};
    keybind goto_next_mo{.key = GLFW_KEY_RIGHT};
    keybind goto_prev_mo{.key = GLFW_KEY_LEFT};
    keybind goto_next_mo_alt{.is_mouse = true, .key = GLFW_MOUSE_BUTTON_5};
    keybind goto_prev_mo_alt{.is_mouse = true, .key = GLFW_MOUSE_BUTTON_4};
    keybind open_config{.key = GLFW_KEY_T};
    keybind rerender{.key = GLFW_KEY_R};
    keybind toggle_render_mo{.key = GLFW_KEY_Z};
    keybind toggle_render_atoms{.key = GLFW_KEY_X};
    keybind toggle_render_box{.key = GLFW_KEY_C};
    keybind save_solution{.key = GLFW_KEY_S, .mods = GLFW_MOD_CONTROL};

    void save(const std::string& path) const;
    void load(const std::string& path);
};

auto on_keybind(const keybind& keybind, int key, int mods, const std::function<void()>& handler, bool mouse) -> bool;

inline static constexpr auto GLFW_MASK_RELEVENT = GLFW_MOD_ALT | GLFW_MOD_CONTROL | GLFW_MOD_SHIFT | GLFW_MOD_SUPER;
