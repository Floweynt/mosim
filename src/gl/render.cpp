#include "gl/render.h"
#include <fmt/core.h>
#include <glm/gtc/type_ptr.hpp>

#define make_shader_set_impl(func, type)                                                                                                             \
    auto gl::shader_manager::func(const std::string& name, type value)->gl::shader_manager&                                                          \
    {                                                                                                                                                \
        expect((bool)current_shader, "current shader must be non-null to set uniform");                                                              \
        current_shader->func(name, value);                                                                                                           \
        return *this;                                                                                                                                \
    }

make_shader_set_impl(set_int, std::int32_t);
make_shader_set_impl(set_vec_int, i32vec2);
make_shader_set_impl(set_vec_int, i32vec3);
make_shader_set_impl(set_vec_int, i32vec4);
make_shader_set_impl(set_float, float);
make_shader_set_impl(set_vec_float, glm::vec2);
make_shader_set_impl(set_vec_float, glm::vec3);
make_shader_set_impl(set_vec_float, glm::vec4);
make_shader_set_impl(set_uint, std::uint32_t);
make_shader_set_impl(set_vec_uint, u32vec2);
make_shader_set_impl(set_vec_uint, u32vec3);
make_shader_set_impl(set_vec_uint, u32vec4);
make_shader_set_impl(set_matrix, const glm::mat2&);
make_shader_set_impl(set_matrix, const glm::mat3&);
make_shader_set_impl(set_matrix, const glm::mat4&);

auto gl::shader_manager::add_shader(const std::string& name, std::string vertex_path, std::string fragment_path, std::string geom_path)
    -> gl::shader_manager&
{
    expect_false(shaders.contains(name), fmt::format("Failed to add shader '{}', it already exists", name));
    shaders.emplace(name, std::make_shared<shader>(std::move(vertex_path), std::move(fragment_path), std::move(geom_path)));
    return *this;
}

auto gl::shader_manager::delete_shader(const std::string& name) -> gl::shader_manager&
{
    expect(shaders.contains(name), fmt::format("Failed to delete nonexistent shader '{}'", name));
    shaders.erase(name);
    return *this;
}

auto gl::shader_manager::use_shader(const std::string& name) -> gl::shader_manager&
{
    expect(shaders.contains(name), fmt::format("Failed to use nonexistent shader '{}'", name));
    current_shader = shaders.at(name);
    current_shader->use();
    return *this;
}

auto gl::shader_manager::push_shader(const std::string& new_shader) -> gl::shader_manager&
{
    expect(shaders.contains(new_shader), fmt::format("Failed to switch to nonexistent shader '{}'", new_shader));
    shader_stack.push(std::move(current_shader));
    current_shader = shaders.at(new_shader);
    current_shader->use();
    return *this;
}

auto gl::shader_manager::pop_shader() -> gl::shader_manager&
{
    expect_false(shader_stack.empty(), "Failed to pop shader stack, it is empty");
    current_shader = std::move(shader_stack.top());
    shader_stack.pop();
    current_shader->use();
    return *this;
}

auto gl::shader_manager::compile_shader(const std::string& name) -> gl::shader_manager&
{
    expect(shaders.contains(name), fmt::format("Failed to compile nonexistent shader '{}'", name));
    auto shader = shaders.at(name);
    shader->compile();
    return *this;
}

auto gl::shader_manager::compile_shader(const std::string& name, const std::function<void(std::string)>& on_error) -> gl::shader_manager&
{
    expect(shaders.contains(name), fmt::format("Failed to compile nonexistent shader '{}'", name));
    shaders.at(name)->compile(on_error);
    return *this;
}

auto gl::shader_manager::compile_shader(const std::string& name, const std::function<void(std::string)>& on_compile_error,
                                        const std::function<void(std::string)>& on_link_error) -> gl::shader_manager&
{
    expect(shaders.contains(name), fmt::format("Failed to compile nonexistent shader '{}'", name));
    shaders.at(name)->compile(on_compile_error, on_link_error);
    return *this;
}

