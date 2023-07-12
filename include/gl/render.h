#pragma once

#include "matrix.h"
#include "shader.h"
#include "util.h"
#include "vertex_buffer.h"
#include <cstdint>
#include <glm/glm.hpp>
#include <memory>
#include <stack>
#include <stdexcept>

namespace gl
{
    class shader_manager
    {
        std::unordered_map<std::string, std::shared_ptr<shader>> shaders;
        std::stack<std::shared_ptr<shader>> shader_stack;
        std::shared_ptr<shader> current_shader;

    public:
        auto has_shader(const std::string& name) -> bool { return shaders.contains(name); }
        auto add_shader(const std::string& name, std::string vertex_path, std::string fragment_path, std::string geom_path = "") -> shader_manager&;

        auto compile_shader(const std::string& name) -> shader_manager&;
        auto compile_shader(const std::string& name, const std::function<void(std::string)>& on_error) -> shader_manager&;
        auto compile_shader(const std::string& name, const std::function<void(std::string)>& on_compile_error,
                                       const std::function<void(std::string)>& on_link_error) -> shader_manager&;

        auto delete_shader(const std::string& name) -> shader_manager&;
        auto use_shader(const std::string& name) -> shader_manager&;

        auto push_shader(const std::string& new_shader) -> shader_manager&;
        auto pop_shader() -> shader_manager&;

        auto get_current_shader() const -> std::shared_ptr<shader> { return current_shader; }

        auto set_int(const std::string& name, std::int32_t value) -> gl::shader_manager&;
        auto set_vec_int(const std::string& name, vec<2, std::int32_t> value) -> gl::shader_manager&;
        auto set_vec_int(const std::string& name, vec<3, std::int32_t> value) -> gl::shader_manager&;
        auto set_vec_int(const std::string& name, vec<4, std::int32_t> value) -> gl::shader_manager&;

        auto set_float(const std::string& name, float value) -> gl::shader_manager&;
        auto set_vec_float(const std::string& name, vec2 value) -> gl::shader_manager&;
        auto set_vec_float(const std::string& name, vec3 value) -> gl::shader_manager&;
        auto set_vec_float(const std::string& name, vec4 value) -> gl::shader_manager&;

        auto set_uint(const std::string& name, std::uint32_t value) -> gl::shader_manager&;
        auto set_vec_uint(const std::string& name, vec<2, std::uint32_t> value) -> gl::shader_manager&;
        auto set_vec_uint(const std::string& name, vec<3, std::uint32_t> value) -> gl::shader_manager&;
        auto set_vec_uint(const std::string& name, vec<4, std::uint32_t> value) -> gl::shader_manager&;

        auto set_matrix(const std::string& name, const mat2& matrix) -> gl::shader_manager&;
        auto set_matrix(const std::string& name, const mat3& matrix) -> gl::shader_manager&;
        auto set_matrix(const std::string& name, const mat4& matrix) -> gl::shader_manager&;
    };

    class render_manager
    {
        shader_manager shader_man;
        matrix_stack model_stack;
        matrix_stack projection_stack;
        matrix_stack view_stack;

        uvec2 screen_size;
    public:
        constexpr void set_screen_size(uvec2 s) { screen_size = s; }
        constexpr auto get_screen_size() -> uvec2 { return screen_size; }

        [[nodiscard]] auto get_shader_manager() -> shader_manager& { return shader_man; }
        [[nodiscard]] auto model() -> matrix_stack& { return model_stack; }
        [[nodiscard]] auto projection() -> matrix_stack& { return projection_stack; }
        [[nodiscard]] auto view() -> matrix_stack& { return view_stack; }

        template <vertex_buffer_cfg T>
        void render_vertex_buffer(const vertex_buffer<T>& buffer)
        {
            expect(buffer.baked(), "buffers must be baked before rendering");
            expect((bool)get_shader_manager().get_current_shader(), "need a shader before rendering");

            shader_man.set_matrix("model", model_stack.curr());
            shader_man.set_matrix("projection", projection_stack.curr());
            shader_man.set_matrix("view", view_stack.curr());
            buffer.render();
        }

        void render_generic(auto functor)
        {
            shader_man.set_matrix("model", model_stack.curr());
            shader_man.set_matrix("projection", projection_stack.curr());
            shader_man.set_matrix("view", view_stack.curr());
            functor();
        }
    };
} // namespace gl
