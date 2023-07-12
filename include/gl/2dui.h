#pragma once
#include "gl/event.h"
#include "gl/font.h"
#include "gl/glstate.h"
#include "gl/render.h"
#include "gl/vertex_buffer.h"
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace gl
{
    // 2d UI component
    class ui_widget
    {
        using buffer_2d = vertex_buffer<vertex_buffer_cfg{.enable_color = true, .primitive = GL_TRIANGLES}>;

        bool dirty{true};
        buffer_2d buffer;
        glm::vec2 ui_pos{};
        glm::vec2 size{};
        std::string name;
        std::string type;
        std::vector<std::tuple<std::string, glm::vec2, float, glm::vec3, int, int>> strings_to_draw;
        glm::uvec2 prev_win_size{};

        std::vector<std::shared_ptr<ui_widget>> children;

        void update_dirty_flag(glm::uvec2 win_size);
        void do_update(glm::uvec2 win_size);

    protected:
        virtual void update(glm::uvec2 win_size) = 0;
        virtual void render_hook(render_manager& renderer) {}

        void draw_text(glm::vec2 pos, const std::string& text, glm::vec3 color, float scale = 1, int horiz = TEXT_LEFT, int vert = TEXT_TOP);
        void draw_rect(glm::vec2 start, glm::vec2 end, glm::vec3 color);
        void draw_rect_border(glm::vec2 start, glm::vec2 end, glm::vec3 color, float thickness = 1);
        void draw_triangle(glm::vec2 vert1, glm::vec2 vert2, glm::vec2 vert3, glm::vec3 color);
        void draw_line(glm::vec2 from, glm::vec2 to, glm::vec3 color, float thickness = 1);

        static constexpr auto in_bb_corners(glm::vec2 from, glm::vec2 to, glm::vec2 target) -> bool
        {
            swap_smaller_first(from.x, to.x);
            swap_smaller_first(from.y, to.y);
            return from.x < target.x && target.x < to.x && from.y < target.y && target.y < to.y;
        }

        static constexpr auto in_bb_size(glm::vec2 from, glm::vec2 size, glm::vec2 target) -> bool
        {
            return in_bb_corners(from, from + size, target);
        }

        constexpr void mark_dirty() { dirty = true; }

        void add_child(std::shared_ptr<ui_widget> child);
        void set_child(std::shared_ptr<ui_widget> child, size_t index);
        void remove_child(const std::shared_ptr<ui_widget>& child);
        void update_children(glm::uvec2 win_size);

        [[nodiscard]] constexpr auto get_children() const -> const auto& { return children; }

    public:
        ui_widget(std::string type, std::string name) : name(std::move(name)), type(std::move(type)) {}

        void render(render_manager& renderer);
        [[nodiscard]] constexpr auto pos() const -> const auto& { return ui_pos; }
        [[nodiscard]] constexpr auto get_size() const -> const auto& { return size; }

        constexpr void set_size(glm::vec2 new_size)
        {
            if (new_size != size)
            {
                mark_dirty();
                size = new_size;
            }
        }

        constexpr void move(glm::vec2 new_pos)
        {
            if (ui_pos != new_pos)
            {
                mark_dirty();
                ui_pos = new_pos;
            }
        }

        constexpr void move(float x, float y) { move({x, y}); }

        virtual ~ui_widget() = default;
    };

    using widget_ref = std::shared_ptr<ui_widget>;

    class ui_button : public ui_widget, public event_handler
    {
        std::string label;
        glm::vec2 text_offset;
        glm::vec3 bg_color;
        glm::vec3 hover_color;
        glm::vec3 text_color;
        std::function<bool(int key, int mods)> on_press;
        float text_scale;

        bool is_hovering{};

        constexpr void set_hovering(bool value)
        {
            if (value != is_hovering)
            {
                is_hovering = value;
                mark_dirty();
            }
        }

    protected:
        void update(glm::uvec2 /*win_size*/) override;

    public:
        ui_button(std::string name, std::string label, glm::vec2 size, glm::vec2 text_offset, glm::vec3 bg_color, glm::vec3 hover_color,
                  glm::vec3 text_color, std::function<bool(int key, int mods)> on_press, float text_scale);

        auto handle_mouse_press(int key, int mods) -> bool override;
        auto handle_mouse(float x, float y) -> bool override;
    };

    class vertical_stack : public ui_widget, public event_handler
    {
        float margin;

    protected:
        void update(glm::uvec2 win_size) override;
        void render_hook(render_manager& renderer) override;

    public:
        vertical_stack(std::string name, float margin) : ui_widget("vertical_stack", std::move(name)), margin(margin) {}
        using ui_widget::add_child;
        using ui_widget::remove_child;

        auto handle_mouse(float x, float y) -> bool override;
        auto handle_mouse_press(int key, int mods) -> bool override;
        auto handle_mouse_release(int key, int mods) -> bool override;
        auto handle_key_press(int key, int scancode, int mods) -> bool override;
        auto handle_key_release(int key, int scancode, int mods) -> bool override;
        ~vertical_stack() override = default;
    };

    class ui_glue : public ui_widget, public event_handler
    {
    public:
        enum glue_target
        {
            GLUE_TOP_LEFT = 0,
            GLUE_TOP,
            GLUE_TOP_RIGHT,
            GLUE_RIGHT,
            GLUE_BOTTOM_RIGHT,
            GLUE_BOTTOM,
            GLUE_BOTTOM_LEFT,
            GLUE_LEFT
        };

    private:
        glue_target target;
        float pos;

    protected:
        void update(glm::uvec2 win_size) override;
        void render_hook(render_manager& renderer) override;

    public:
        ui_glue(std::string name, glue_target target, widget_ref child);
        ui_glue(std::string name, glue_target target, float pos, widget_ref child);

        auto handle_mouse(float x, float y) -> bool override;
        auto handle_mouse_press(int key, int mods) -> bool override;
        auto handle_mouse_release(int key, int mods) -> bool override;
        auto handle_key_press(int key, int scancode, int mods) -> bool override;
        auto handle_key_release(int key, int scancode, int mods) -> bool override;

        ~ui_glue() override = default;
    };

    class ui_padding : public ui_widget, public event_handler
    {
    private:
        float top;
        float bottom;
        float left;
        float right;

    protected:
        void update(glm::uvec2 win_size) override;
        void render_hook(render_manager& renderer) override;

    public:
        ui_padding(std::string name, float top, float bottom, float left, float right, widget_ref child);

        auto handle_mouse(float x, float y) -> bool override;
        auto handle_mouse_press(int key, int mods) -> bool override;
        auto handle_mouse_release(int key, int mods) -> bool override;
        auto handle_key_press(int key, int scancode, int mods) -> bool override;
        auto handle_key_release(int key, int scancode, int mods) -> bool override;

        ~ui_padding() override = default;
    };

    class ui_label : public ui_widget
    {
    private:
        std::string text;
        glm::vec2 text_offset;
        glm::vec3 bg_color;
        glm::vec3 text_color;
        float text_scale;

    protected:
        void update(glm::uvec2 win_size) override;

    public:
        ui_label(std::string name, std::string text, glm::vec3 bg_color, glm::vec3 text_color, glm::vec2 text_offset, float text_scale);
        ~ui_label() override = default;
    };
} // namespace gl
