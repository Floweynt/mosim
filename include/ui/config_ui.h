#pragma once

// Label     [keybind] alt [keybind]
#include "gl/2dui.h"
#include "gl/2dui_builder.h"
#include "gl/event.h"
#include "gl/font.h"
#include "keybind.h"
#include "platform_dep.h"
#include "ui_config.h"
#include <GLFW/glfw3.h>
#include <fmt/core.h>
#include <memory>
#include <string>

class keybind_ui : public gl::ui_widget, public gl::event_handler
{
    std::string label;

    keybind* main;
    keybind* alt;
    keybind* keybind_target{};
    keybind* hovering_target{};
    std::function<void()> on_update;

    void set_keybind_target(keybind* val);
    void set_hovering_target(keybind* val);

    inline static constexpr auto MAIN_KB_POS_1 = glm::vec2(MAIN_KB_OFFSET, (KB_UI_SIZE.y - KB_DISPLAY_SIZE.y) / 2);
    inline static constexpr auto MAIN_KB_POS_2 = glm::vec2(MAIN_KB_OFFSET + KB_DISPLAY_SIZE.x, (KB_UI_SIZE.y + KB_DISPLAY_SIZE.y) / 2);
    inline static constexpr auto ALT_KB_POS_1 = glm::vec2(ALT_KB_OFFSET, (KB_UI_SIZE.y - KB_DISPLAY_SIZE.y) / 2);
    inline static constexpr auto ALT_KB_POS_2 = glm::vec2(ALT_KB_OFFSET + KB_DISPLAY_SIZE.x, (KB_UI_SIZE.y + KB_DISPLAY_SIZE.y) / 2);

protected:
    void update(glm::uvec2 win_size) override;

public:
    auto handle_mouse(float x, float y) -> bool override;
    auto handle_mouse_press(int key, int mods) -> bool override;
    auto handle_key_press(int key, int scancode, int mods) -> bool override;

    keybind_ui(std::string name, std::string label, std::function<void()> on_update, keybind* main, keybind* alt = nullptr);
};

class config_ui : public gl::ui_widget, public gl::event_handler
{
    gl::widget_ref keybind_stack;
    bool is_active{};

protected:
    void update(glm::uvec2 win_size) override;
    void render_hook(gl::render_manager& renderer) override;
    static void on_update();

public:
    config_ui(std::string name);

    auto handle_mouse(float x, float y) -> bool override;
    auto handle_mouse_press(int key, int mods) -> bool override;
    auto handle_mouse_release(int key, int mods) -> bool override;
    auto handle_key_press(int key, int scancode, int mods) -> bool override;
    auto handle_key_release(int key, int scancode, int mods) -> bool override;
    void set_active(bool flag);
};
