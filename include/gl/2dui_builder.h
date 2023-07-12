#pragma once
#include "gl/2dui.h"

#define UI_BUILDER_MAKE_BUILD_METHOD(var_type, method_name, var_name)                                                                                \
    constexpr auto method_name(var_type val)->auto&                                                                                                  \
    {                                                                                                                                                \
        var_name = val;                                                                                                                              \
        return *this;                                                                                                                                \
    }

namespace gl
{
    class button_builder
    {
        std::string name;
        std::string label;
        glm::vec2 _size;
        glm::vec2 _text_offset;
        glm::vec3 _bg_color;
        glm::vec3 _hover_color;
        glm::vec3 _text_color;
        std::function<bool(int, int)> _on_press;
        float _text_scale{1};

    public:
        button_builder(std::string name, std::string label)
            : name(std::move(name)), label(std::move(label)), _size(10, 10), _text_offset(0, 0), _bg_color(0, 0, 0), _hover_color(0, 0, 0),
              _text_color(1, 1, 1)
        {
        }

        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec2, size, _size);
        constexpr auto size(float width, float height) -> auto& { return size({width, height}); }

        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec2, text_offset, _text_offset);
        constexpr auto text_offset(float x_off, float y_off) -> auto& { return text_offset({x_off, y_off}); }

        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec3, bg_color, _bg_color);
        constexpr auto bg_color(float r, float g, float b) -> auto& { return bg_color({r, g, b}); }

        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec3, hover_color, _hover_color);
        constexpr auto hover_color(float r, float g, float b) -> auto& { return hover_color({r, g, b}); }

        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec3, text_color, _text_color);
        constexpr auto text_color(float r, float g, float b) -> auto& { return text_color({r, g, b}); }

        UI_BUILDER_MAKE_BUILD_METHOD(float, text_scale, _text_scale);

        auto on_press(std::function<bool(int, int)> val) -> button_builder&
        {
            _on_press = std::move(val);
            return *this;
        }

        auto build() -> widget_ref
        {
            return std::make_shared<ui_button>(std::move(name), std::move(label), _size, _text_offset, _bg_color, _hover_color, _text_color,
                                               std::move(_on_press), _text_scale);
        }
    };

    class vertical_stack_builder
    {
        std::string name;
        float _margin{};
        std::vector<widget_ref> children;

    public:
        vertical_stack_builder(std::string name) : name(std::move(name)) {}
        UI_BUILDER_MAKE_BUILD_METHOD(float, margin, _margin);

        auto add_child(widget_ref ref) -> vertical_stack_builder&
        {
            children.emplace_back(std::move(ref));
            return *this;
        }

        auto build() -> widget_ref
        {
            auto res = std::make_shared<vertical_stack>(std::move(name), _margin);
            for (auto& child : children)
            {
                res->add_child(std::move(child));
            }
            return res;
        }
    };

    inline auto build_glue(std::string name, ui_glue::glue_target target, widget_ref child) -> widget_ref
    {
        return std::make_shared<ui_glue>(std::move(name), target, std::move(child));
    }

    inline auto build_glue(std::string name, float pos, ui_glue::glue_target target, widget_ref child) -> widget_ref
    {
        return std::make_shared<ui_glue>(std::move(name), target, pos, std::move(child));
    }

    class padding_builder
    {
        std::string name;
        float _top{};
        float _bottom{};
        float _left{};
        float _right{};
        widget_ref child;

    public:
        padding_builder(std::string name, widget_ref child) : name(std::move(name)), child(std::move(child)) {}
        UI_BUILDER_MAKE_BUILD_METHOD(float, top, _top);
        UI_BUILDER_MAKE_BUILD_METHOD(float, bottom, _bottom);
        UI_BUILDER_MAKE_BUILD_METHOD(float, left, _left);
        UI_BUILDER_MAKE_BUILD_METHOD(float, right, _right);

        constexpr auto pad(float val) -> auto&
        {
            _top = _left = _bottom = _right = val;
            return *this;
        }

        auto build() -> widget_ref
        {
            auto res = std::make_shared<ui_padding>(std::move(name), _top, _bottom, _left, _right, std::move(child));
            return res;
        }
    };

    class label_builder
    {
        std::string name;
        std::string text;
        glm::vec2 _text_offset{0, 0};
        glm::vec3 _bg_color{0, 0, 0};
        glm::vec3 _text_color{1, 1, 1};
        float _text_scale{1};
        glm::vec2 _size{10, 10};

    public:
        label_builder(std::string name, std::string text) : name(std::move(name)), text(std::move(text)) {}
        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec2, size, _size);
        constexpr auto size(float x, float y) -> auto& { return size({x, y}); }
        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec2, text_offset, _text_offset);
        constexpr auto text_offset(float x, float y) -> auto& { return text_offset({x, y}); }
        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec3, bg_color, _bg_color);
        constexpr auto bg_color(float r, float g, float b) -> auto& { return bg_color({r, g, b}); }
        UI_BUILDER_MAKE_BUILD_METHOD(glm::vec3, text_color, _text_color);
        constexpr auto text_color(float r, float g, float b) -> auto& { return text_color({r, g, b}); }
        UI_BUILDER_MAKE_BUILD_METHOD(float, text_scale, _text_scale);

        auto build() -> widget_ref
        {
            auto res = std::make_shared<ui_label>(std::move(name), std::move(text), _bg_color, _text_color, _text_offset, _text_scale);
            res->set_size(_size);
            return res;
        }
    };
} // namespace gl

#undef UI_BUILDER_MAKE_BUILD_METHOD
