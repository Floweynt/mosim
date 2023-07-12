#include "gl/event.h"
using namespace gl;

auto compound_event_handler::handle_mouse(float x, float y) -> bool
{
    for (const auto& handler : handlers)
    {
        if (handler->handle_mouse(x, y))
        {
            return true;
        }
    }
    return mouse_handler ? mouse_handler(x, y) : false;
}

auto compound_event_handler::handle_mouse_press(int key, int mods) -> bool
{
    for (const auto& handler : handlers)
    {
        if (handler->handle_mouse_press(key, mods))
        {
            return true;
        }
    }
    return mouse_press_handler ? mouse_press_handler(key, mods) : false;
}

auto compound_event_handler::handle_mouse_release(int key, int mods) -> bool
{
    for (const auto& handler : handlers)
    {
        if (handler->handle_mouse_release(key, mods))
        {
            return true;
        }
    }
    return mouse_release_handler ? mouse_release_handler(key, mods) : false;
}

auto compound_event_handler::handle_key_press(int key, int scancode, int mods) -> bool
{
    for (const auto& handler : handlers)
    {
        if (handler->handle_key_press(key, scancode, mods))
        {
            return true;
        }
    }
    return key_press_handler ? key_press_handler(key, scancode, mods) : false;
}

auto compound_event_handler::handle_key_release(int key, int scancode, int mods) -> bool
{
    for (const auto& handler : handlers)
    {
        if (handler->handle_key_release(key, scancode, mods))
        {
            return true;
        }
    }
    return key_release_handler ? key_release_handler(key, scancode, mods) : false;
}

