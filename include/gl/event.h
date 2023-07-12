#pragma once

#include "util.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>
namespace gl
{
    class event_handler
    {
    public:
        virtual auto handle_mouse(float /*x*/, float /*y*/) -> bool { return false; }
        virtual auto handle_mouse_press(int /*key*/, int /*mods*/) -> bool { return false; }
        virtual auto handle_mouse_release(int /*key*/, int /*mods*/) -> bool { return false; }
        virtual auto handle_key_press(int /*key*/, int /*scancode*/, int /*mods*/) -> bool { return false; }
        virtual auto handle_key_release(int /*key*/, int /*scancode*/, int /*mods*/) -> bool { return false; }
        virtual ~event_handler() = default;
    };

    class compound_event_handler : public event_handler
    {
        std::vector<std::shared_ptr<event_handler>> handlers;
        std::function<bool(float, float)> mouse_handler;
        std::function<bool(int, int)> mouse_press_handler;
        std::function<bool(int, int)> mouse_release_handler;
        std::function<bool(int, int, int)> key_press_handler;
        std::function<bool(int, int, int)> key_release_handler;

    protected:
        template <typename T>
        void add_child(std::shared_ptr<T> handler)
        {
            std::shared_ptr<event_handler> handler_real = std::dynamic_pointer_cast<event_handler>(handler);
            expect((bool)handler_real, "cannot add child handler that is null");
            handlers.emplace_back(std::move(handler_real));
        }

        void set_mouse_handler(std::function<bool(float, float)> handler) { mouse_handler = std::move(handler); }
        void set_mouse_press_handler(std::function<bool(int, int)> handler) { mouse_press_handler = std::move(handler); }
        void set_mouse_release_handler(std::function<bool(int, int)> handler) { mouse_release_handler = std::move(handler); }
        void set_key_press_handler(std::function<bool(int, int, int)> handler) { key_press_handler = std::move(handler); }
        void set_key_release_handler(std::function<bool(int, int, int)> handler) { key_release_handler = std::move(handler); }

    public:
        auto handle_mouse(float x, float y) -> bool override;
        auto handle_mouse_press(int key, int mods) -> bool override;
        auto handle_mouse_release(int key, int mods) -> bool override;
        auto handle_key_press(int key, int scancode, int mods) -> bool override;
        auto handle_key_release(int key, int scancode, int mods) -> bool override;
    };

    template <int type>
    class mouse_drag_handler : public event_handler
    {
        bool is_pressed = false;

        float prev_mouse_x = NAN;
        float prev_mouse_y{};

        std::function<bool(float, float)> drag_delta_callback;

    public:
        auto handle_mouse(float mouse_x, float mouse_y) -> bool override
        {
            float delta_x = mouse_x - prev_mouse_x;
            float delta_y = mouse_y - prev_mouse_y;
            bool ret = false;
            if (delta_x != NAN)
            {
                if (is_pressed)
                {
                    ret = drag_delta_callback ? drag_delta_callback(delta_x, delta_y) : false;
                }
            }

            prev_mouse_x = mouse_x;
            prev_mouse_y = mouse_y;
            return ret;
        }

        auto handle_mouse_press(int key, int /*unused*/) -> bool override
        {
            if (key == type)
            {
                is_pressed = true;
                return true;
            }
            return false;
        }

        auto handle_mouse_release(int key, int /*unused*/) -> bool override
        {
            if (key == type)
            {
                is_pressed = false;
                return true;
            }
            return false;
        }

        void set_handler(std::function<bool(float, float)> handler) { drag_delta_callback = std::move(handler); }
    };
} // namespace gl
