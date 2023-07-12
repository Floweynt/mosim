#include "gl/2dui.h"
#include "gl/font.h"

using namespace gl;

void ui_widget::update_dirty_flag(glm::uvec2 win_size)
{
    if (dirty)
    {
        return;
    }

    if (win_size != prev_win_size)
    {
        mark_dirty();
        return;
    }

    for (const auto& child : children)
    {
        child->update_dirty_flag(win_size);
    }

    for (const auto& child : children)
    {
        if (child->dirty)
        {
            mark_dirty();
            return;
        }
    }
}

void ui_widget::do_update(glm::uvec2 win_size)
{
    update_dirty_flag(win_size);
    if (!dirty)
    {
        return;
    }

    prev_win_size = win_size;

    while (dirty)
    {
        update_children(win_size);
        strings_to_draw.clear();
        buffer.reset();
#ifdef DEBUG
        std::cout << "update_tracker: " << type << "::" << name << "\n";
#endif
        dirty = false;
        update(win_size);
        update_dirty_flag(win_size);
    }

    buffer.bake();
}

void ui_widget::add_child(std::shared_ptr<ui_widget> child)
{
    expect(child.get() != this, "cannot add current widget as it's own child");
    children.emplace_back(std::move(child));
    mark_dirty();
}

void ui_widget::set_child(std::shared_ptr<ui_widget> child, size_t index)
{
    expect(child.get() != this, "cannot add current widget as it's own child");
    children[index] = std::move(child);
    mark_dirty();
}

void ui_widget::remove_child(const std::shared_ptr<ui_widget>& child) { std::erase(children, child); }

void ui_widget::update_children(glm::uvec2 win_size)
{
    for (const auto& child : children)
    {
        child->do_update(win_size);
    }
}

void ui_widget::draw_text(glm::vec2 pos, const std::string& text, glm::vec3 color, float scale, int horiz, int vert)
{
    strings_to_draw.emplace_back(text, pos + ui_pos, scale, color, horiz, vert);
}

void ui_widget::draw_rect(glm::vec2 start, glm::vec2 end, glm::vec3 color)
{
    // s, s   e, s
    // +------+
    // |      |
    // +------+
    // s, e   e e
    start += ui_pos;
    end += ui_pos;
    buffer.vert(buffer_2d::builder().vert(start.x, start.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(start.x, end.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(end.x, start.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(end.x, end.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(start.x, end.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(end.x, start.y, 0).colorf(color.r, color.g, color.b, 1).end());
}

void ui_widget::draw_triangle(glm::vec2 vert1, glm::vec2 vert2, glm::vec2 vert3, glm::vec3 color)
{
    vert1 += ui_pos;
    vert2 += ui_pos;
    vert3 += ui_pos;

    buffer.vert(buffer_2d::builder().vert(vert1.x, vert1.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(vert2.x, vert2.y, 0).colorf(color.r, color.g, color.b, 1).end());
    buffer.vert(buffer_2d::builder().vert(vert3.x, vert3.y, 0).colorf(color.r, color.g, color.b, 1).end());
}

void ui_widget::draw_line(glm::vec2 from, glm::vec2 to, glm::vec3 color, float thickness)
{
    from += ui_pos;
    to += ui_pos;
    glm::vec2 disp = to - from;
    glm::vec2 prep_vec = glm::normalize(glm::vec2{disp.y, disp.x});
    thickness /= 2;

    draw_rect(from + prep_vec * thickness, to - prep_vec * thickness, color);
}

void ui_widget::render(render_manager& renderer)
{
    do_update(renderer.get_screen_size());
    if (!renderer.get_shader_manager().has_shader("2d_shader"))
    {
        renderer.get_shader_manager()
            .add_shader("2d_shader", "assets/shaders/2d_vertex.glsl", "assets/shaders/2d_fragment.glsl")
            .compile_shader("2d_shader");
    }

    renderer.get_shader_manager()
        .push_shader("2d_shader")
        .set_matrix("projection", glm::ortho(0.0f, (float)renderer.get_screen_size().x, (float)renderer.get_screen_size().y, 0.0f, -1.0f, 1.0f));

    glDisable(GL_DEPTH_TEST);
    gl::gl_check();
    buffer.get_vao().draw(GL_TRIANGLES, 0, buffer.vert_count());
    render_hook(renderer);
    glEnable(GL_DEPTH_TEST);
    gl::gl_check();

    renderer.get_shader_manager().pop_shader();

    for (const auto& [string, pos, scale, color, horiz, vert] : strings_to_draw)
    {
        render_text(string, pos, scale, glm::vec4(color, 1), horiz, vert);
    }
}
void ui_button::update(glm::uvec2 /*win_size*/)
{
    draw_rect({0, 0}, get_size(), is_hovering ? hover_color : bg_color);
    draw_text(get_size() / 2.f + text_offset, label, text_color, text_scale, TEXT_CENTER, TEXT_CENTER);
}

ui_button::ui_button(std::string name, std::string label, glm::vec2 size, glm::vec2 text_offset, glm::vec3 bg_color, glm::vec3 hover_color,
                     glm::vec3 text_color, std::function<bool(int key, int mods)> on_press, float text_scale)
    : ui_widget("button", std::move(name)), label(std::move(label)), text_offset(text_offset), bg_color(bg_color), hover_color(hover_color),
      text_color(text_color), on_press(std::move(on_press)), text_scale(text_scale)
{
    set_size(size);
}

auto ui_button::handle_mouse_press(int key, int mods) -> bool
{
    if (!is_hovering)
    {
        return false;
    }

    return (bool)on_press ? on_press(key, mods) : false;
}

auto ui_button::handle_mouse(float x, float y) -> bool
{
    glm::vec2 mouse_pos(x, y);
    if (in_bb_size(pos(), get_size(), mouse_pos))
    {
        set_hovering(true);
        return true;
    }
    set_hovering(false);
    return false;
}

void vertical_stack::update(glm::uvec2 win_size)
{
    float offset = 0;
    float max_x = 0;
    for (const auto& child : get_children())
    {
        child->move(pos().x, pos().y + offset);
        max_x = std::max(child->get_size().x, max_x);
        offset += child->get_size().y + margin;
    }

    // we need to subtract the margin since we dont need that on the last element
    set_size({max_x, offset - margin});
}

void vertical_stack::render_hook(render_manager& renderer)
{
    for (const auto& child : get_children())
    {
        child->render(renderer);
    }
}

auto vertical_stack::handle_mouse(float x, float y) -> bool
{
    for (const auto& child : get_children())
    {
        auto ptr = dynamic_pointer_cast<event_handler>(child);
        if (ptr)
        {
            if (ptr->handle_mouse(x, y))
            {
                return true;
            }
        }
    }
    return false;
}

auto vertical_stack::handle_mouse_press(int key, int mods) -> bool
{
    for (const auto& child : get_children())
    {
        auto ptr = dynamic_pointer_cast<event_handler>(child);
        if (ptr)
        {
            if (ptr->handle_mouse_press(key, mods))
            {
                return true;
            }
        }
    }
    return false;
}

auto vertical_stack::handle_mouse_release(int key, int mods) -> bool
{
    for (const auto& child : get_children())
    {
        auto ptr = dynamic_pointer_cast<event_handler>(child);
        if (ptr)
        {
            if (ptr->handle_mouse_release(key, mods))
            {
                return true;
            }
        }
    }
    return false;
}

auto vertical_stack::handle_key_press(int key, int scancode, int mods) -> bool
{
    for (const auto& child : get_children())
    {
        auto ptr = dynamic_pointer_cast<event_handler>(child);
        if (ptr)
        {
            if (ptr->handle_key_press(key, scancode, mods))
            {
                return true;
            }
        }
    }
    return false;
}

auto vertical_stack::handle_key_release(int key, int scancode, int mods) -> bool
{
    for (const auto& child : get_children())
    {
        auto ptr = dynamic_pointer_cast<event_handler>(child);
        if (ptr)
        {
            if (ptr->handle_key_release(key, scancode, mods))
            {
                return true;
            }
        }
    }
    return false;
}

ui_glue::ui_glue(std::string name, glue_target target, widget_ref child) : ui_widget("glue", std::move(name)), target(target), pos(0)
{
    expect_false((target & 1) != 0, "must use corners if using the constructor without position");
    add_child(std::move(child));
}

ui_glue::ui_glue(std::string name, glue_target target, float pos, widget_ref child) : ui_widget("glue", std::move(name)), target(target), pos(pos)
{
    expect_true((target & 1) != 0, "do not specify pos when constructing with a corner");
    add_child(std::move(child));
}

void ui_glue::update(glm::uvec2 win_size)
{
    update_children(win_size);
    // we must glue the object back to the proper corner/wall
    const auto& child = get_children()[0];
    glm::vec2 new_pos;

    switch (target)
    {
    case GLUE_TOP_LEFT:
        new_pos = {0, 0};
        break;
    case GLUE_TOP:
        new_pos = {pos, 0};
        break;
    case GLUE_TOP_RIGHT:
        new_pos = {win_size.x - child->get_size().x, 0};
        break;
    case GLUE_RIGHT:
        new_pos = {win_size.x - child->get_size().x, pos};
        break;
    case GLUE_BOTTOM_RIGHT:
        new_pos = {win_size.x - child->get_size().x, win_size.y - child->get_size().y};
        break;
    case GLUE_BOTTOM:
        new_pos = {pos, win_size.y - child->get_size().y};
        break;
    case GLUE_BOTTOM_LEFT:
        new_pos = {0, win_size.y - child->get_size().y};
        break;
    case GLUE_LEFT:
        new_pos = {0, pos};
        break;
    }

    child->move(new_pos);
    update_children(win_size);

    move(new_pos);
    set_size(child->get_size());
}

void ui_glue::render_hook(render_manager& renderer) { get_children()[0]->render(renderer); }

auto ui_glue::handle_mouse(float x, float y) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse(x, y);
    }
    return false;
}

auto ui_glue::handle_mouse_press(int key, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse_press(key, mods);
    }
    return false;
}

auto ui_glue::handle_mouse_release(int key, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse_release(key, mods);
    }
    return false;
}

auto ui_glue::handle_key_press(int key, int scancode, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_key_press(key, scancode, mods);
    }
    return false;
}

auto ui_glue::handle_key_release(int key, int scancode, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_key_release(key, scancode, mods);
    }
    return false;
}

ui_padding::ui_padding(std::string name, float top, float bottom, float left, float right, widget_ref child)
    : ui_widget("padding", std::move(name)), top(top), bottom(bottom), left(left), right(right)
{
    add_child(std::move(child));
}

void ui_padding::update(glm::uvec2 win_size)
{
    const auto& child = get_children()[0];
    glm::vec2 new_pos = pos() + glm::vec2(left, top);
    child->move(new_pos);
    set_size(child->get_size() + glm::vec2(left + right, top + bottom));
}

void ui_padding::render_hook(render_manager& renderer) { get_children()[0]->render(renderer); }

auto ui_padding::handle_mouse(float x, float y) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse(x, y);
    }
    return false;
}

auto ui_padding::handle_mouse_press(int key, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse_press(key, mods);
    }
    return false;
}

auto ui_padding::handle_mouse_release(int key, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_mouse_release(key, mods);
    }
    return false;
}

auto ui_padding::handle_key_press(int key, int scancode, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_key_press(key, scancode, mods);
    }
    return false;
}

auto ui_padding::handle_key_release(int key, int scancode, int mods) -> bool
{
    auto handler = std::dynamic_pointer_cast<event_handler>(get_children()[0]);
    if (handler)
    {
        return handler->handle_key_release(key, scancode, mods);
    }
    return false;
}

void ui_label::update(glm::uvec2 win_size)
{
    draw_rect({0, 0}, get_size(), bg_color);
    draw_text(get_size() / 2.f + text_offset, text, text_color, text_scale, TEXT_CENTER, TEXT_CENTER);
}

ui_label::ui_label(std::string name, std::string text, glm::vec3 bg_color, glm::vec3 text_color, glm::vec2 text_offset, float text_scale)
    : ui_widget("label", name), text(std::move(text)), bg_color(bg_color), text_color(text_color), text_offset(text_offset), text_scale(text_scale)
{
}

