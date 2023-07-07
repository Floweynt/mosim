#ifndef __VERTEX_BUFFER_H__
#define __VERTEX_BUFFER_H__
#include <GL/glew.h>

// reorder barrier
#include "glstate.h"
#include "matrix.h"
#include "shader.h"
#include <GL/gl.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <glm/gtx/string_cast.hpp>
#include <glm/matrix.hpp>
#include <glm/vec3.hpp>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace gl
{
    using namespace glm;

    struct vertex_buffer_cfg
    {
        bool enable_normal;
        bool enable_uv;
        bool enable_color;
        GLenum primitive = GL_TRIANGLES;
    };

    template <vertex_buffer_cfg C>
    class vertex_buffer
    {
        struct vertex_with_normal
        {
            vec3 norm;
        };

        struct vertex_with_uv
        {
            vec2 uv;
        };

        struct vertex_with_color
        {
            std::uint32_t color;
        };

        struct empty0
        {
        };

        struct empty1
        {
        };

        struct empty2
        {
        };

    public:
        struct vertex : std::conditional_t<C.enable_normal, vertex_with_normal, empty0>,
                        std::conditional_t<C.enable_uv, vertex_with_uv, empty1>,
                        std::conditional_t<C.enable_color, vertex_with_color, empty2>
        {
            vec3 pos;
        };

        static_assert(sizeof(vertex) == sizeof(vec3) + (C.enable_normal ? sizeof(vec3) : 0) + (C.enable_uv ? sizeof(vec2) : 0) +
                                            (C.enable_color ? sizeof(std::uint32_t) : 0));

    private:
        std::vector<vertex> buffer;
        bool is_baked{};
        vertex_array_object vao;
        buffer_object vbo;
        friend class render_manager;

        void render() const { vao.draw(C.primitive, 0, buffer.size()); }

    public:
        constexpr vertex_buffer() = default;

        class builder;

        [[nodiscard]] constexpr auto baked() const -> bool { return is_baked; }

        void bake()
        {
            if (is_baked)
            {
                return;
            }

            vbo.init(std::span(buffer), GL_STATIC_DRAW);
            vao.bind_vbo(vbo, 0, 0, sizeof(vertex));

            vao.format(0, decltype(vertex::pos)::length(), GL_FLOAT, false, offsetof(vertex, pos));
            vao.enable_attribute(0);
            vao.set_attribute_buffer_bind(0, 0);

            if constexpr (C.enable_normal)
            {
                vao.format(1, decltype(vertex::norm)::length(), GL_FLOAT, false, offsetof(vertex, norm));
                vao.enable_attribute(1);
                vao.set_attribute_buffer_bind(1, 0);
            }

            if constexpr (C.enable_uv)
            {
                vao.format(2, decltype(vertex::uv)::length(), GL_FLOAT, false, offsetof(vertex, uv));
                vao.enable_attribute(2);
                vao.set_attribute_buffer_bind(2, 0);
            }

            if constexpr (C.enable_color)
            {
                vao.int_format(3, 1, GL_UNSIGNED_INT, offsetof(vertex, color));
                vao.enable_attribute(3);
                vao.set_attribute_buffer_bind(3, 0);
            }

            is_baked = true;
        }

        auto vert(vertex vtx) -> vertex_buffer&
        {
            buffer.emplace_back(vtx);
            return *this;
        }
    };

    template <vertex_buffer_cfg C>
    class vertex_buffer<C>::builder
    {
        vertex v;

    public:
        constexpr auto vert(const vec3& vec) -> auto&
        {
            v.pos = vec;
            return *this;
        }

        constexpr auto vert(float x, float y, float z) -> auto&
        {
            v.pos = {x, y, z};
            return *this;
        }

        constexpr auto normal(const vec3& vec) -> auto& requires(C.enable_normal) {
            v.norm = vec;
            return *this;
        }

        constexpr auto normal(float x, float y, float z) -> auto& requires(C.enable_normal) {
            v.norm = {x, y, z};
            return *this;
        }

        constexpr auto color(std::uint8_t r, std::uint8_t g, std::uint8_t b, std::uint8_t a = 255) -> auto& requires(C.enable_color) {
            v.color = ((std::uint32_t)a << 24) | ((std::uint32_t)r << 16) | ((std::uint32_t)g << 8) | ((std::uint32_t)b << 0);
            return *this;
        }

        constexpr auto colorf(float r, float g, float b, float a = 1.f) -> auto& requires(C.enable_color) {
            return color((std::uint8_t)(std::clamp(r, 0.f, 1.f) * 255), (std::uint8_t)(std::clamp(g, 0.f, 1.f) * 255),
                         (std::uint8_t)(std::clamp(b, 0.f, 1.f) * 255), (std::uint8_t)(std::clamp(a, 0.f, 1.f) * 255));
        }

        constexpr auto color(const uint32_t color) -> auto& requires(C.enable_color)

        {
            v.color = color;
            return *this;
        }

        constexpr auto uv(float u, float v) -> auto& requires(C.enable_uv)

        {
            this->v.uv = {u, v};
            return *this;
        }

        constexpr auto uv(const vec2& uv) -> auto& requires(C.enable_uv) {
            this->v.uv = uv;
            return *this;
        }

        constexpr auto end() -> vertex
        {
            return v;
        }
    };

    using vb_all = vertex_buffer<vertex_buffer_cfg{
        .enable_normal = true,
        .enable_uv = true,
        .enable_color = true,
    }>;
} // namespace gl

#endif
