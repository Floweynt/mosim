
#include "gl/shader.h"
#include "gl/glstate.h"
#include "resources.h"
#include "util.h"
#include <fmt/core.h>
#include <fstream>
#include <functional>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <sstream>

inline static constexpr auto DEFAULT_ERROR_LOG_SIZE = 4096;

static auto check_compile_error(unsigned int shader, const std::function<void(std::string)>& on_fail) -> bool
{
    int success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    gl::gl_check();

    if (!success)
    {
        int len = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
        gl::gl_check();
        if (len <= 0)
        {
            len = DEFAULT_ERROR_LOG_SIZE;
        }

        std::string msg;
        msg.resize(len);
        glGetShaderInfoLog(shader, len, nullptr, msg.data());
        gl::gl_check();
        on_fail(std::move(msg));
        return false;
    }

    return true;
}

static auto check_link_error(unsigned int shader, const std::function<void(std::string)>& on_fail) -> bool
{
    int success = 0;
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    gl::gl_check();

    if (!success)
    {
        int len = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
        gl::gl_check();
        if (len < 0)
        {
            len = DEFAULT_ERROR_LOG_SIZE;
        }

        std::string msg;
        msg.resize(len);
        glGetShaderInfoLog(shader, len, nullptr, msg.data());
        gl::gl_check();
        on_fail(std::move(msg));
        return false;
    }

    return true;
}

auto gl::basic_shader::compile() -> bool
{
    return compile([](std::string str) { throw shader_error("failed to compile shader", std::move(str)); },
                   [](std::string str) { throw shader_error("failed to link shader", std::move(str)); });
}

auto gl::basic_shader::compile(err_hdl_t on_error) -> bool { return compile(on_error, on_error); }

auto gl::basic_shader::compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool { return do_compile(on_compile_error, on_link_error); }

void gl::basic_shader::use() const
{
    if (!id)
    {
        throw std::runtime_error("attempted to use shader before it was (successfully) compiled");
    }

    glUseProgram(id);
    gl::gl_check();
}

void gl::basic_shader::dispose()
{
    if (id)
    {
        glDeleteProgram(id);
        gl::gl_check();
    }
}

gl::basic_shader::~basic_shader() { dispose(); }

auto gl::shader::do_compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool
{
    const auto& fragment_data = resource_manager::get_instance().get_resource(fragment_path);
    std::string fragment_str(fragment_data.begin(), fragment_data.end());
    const char* fragment_buf = fragment_str.c_str();

    const auto& vertex_data = resource_manager::get_instance().get_resource(vertex_path);
    std::string vertex_str(vertex_data.begin(), vertex_data.end());
    const char* vertex_buf = vertex_str.c_str();

    unsigned int vertex = 0;
    unsigned int fragment = 0;
    unsigned int geom = 0;

    vertex = glCreateShader(GL_VERTEX_SHADER);
    gl::gl_check();
    fragment = glCreateShader(GL_FRAGMENT_SHADER);
    gl::gl_check();

    if (!geom_path.empty())
    {
        geom = glCreateShader(GL_GEOMETRY_SHADER);
        gl::gl_check();
    }

    glShaderSource(vertex, 1, &vertex_buf, nullptr);
    gl::gl_check();
    glShaderSource(fragment, 1, &fragment_buf, nullptr);
    gl::gl_check();
    if (geom)
    {
        const auto& geom_data = resource_manager::get_instance().get_resource(geom_path);
        std::string geom_str(geom_data.begin(), geom_data.end());
        const char* geom_buf = geom_str.c_str();
        glShaderSource(geom, 1, &geom_buf, nullptr);
        gl::gl_check();
    }

    glCompileShader(vertex);
    gl::gl_check();
    glCompileShader(fragment);
    gl::gl_check();

    if (geom)
    {
        glCompileShader(geom);
        gl::gl_check();
    }

    bool success = true;
    success &= check_compile_error(vertex, on_compile_error);
    success &= check_compile_error(fragment, on_compile_error);
    if (geom)
    {
        success &= check_compile_error(geom, on_compile_error);
    }

    if (success)
    {
        id = glCreateProgram();
        gl::gl_check();
        glAttachShader(id, vertex);
        gl::gl_check();
        glAttachShader(id, fragment);
        gl::gl_check();

        if (geom)
        {
            glAttachShader(id, geom);
            gl::gl_check();
        }

        glLinkProgram(id);
        gl::gl_check();
        success &= check_link_error(id, on_link_error);
    }

    glDeleteShader(vertex);
    gl::gl_check();
    glDeleteShader(fragment);
    gl::gl_check();

    if (geom)
    {
        glDeleteShader(geom);
        gl::gl_check();
    }

    if (!success)
    {
        id = 0;
    }

    return success;
}

auto gl::compute_shader::do_compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool
{
    const auto& shader_data = resource_manager::get_instance().get_resource(path);
    std::string shader_str(shader_data.begin(), shader_data.end());
    const char* shader_buf = shader_str.c_str();

    unsigned int compute = glCreateShader(GL_COMPUTE_SHADER);
    gl::gl_check();
    glShaderSource(compute, 1, &shader_buf, nullptr);
    gl::gl_check();
    glCompileShader(compute);
    gl::gl_check();

    bool success = true;
    success &= check_compile_error(compute, on_compile_error);

    if (success)
    {
        id = glCreateProgram();
        gl::gl_check();
        glAttachShader(id, compute);
        gl::gl_check();
        glLinkProgram(id);
        gl::gl_check();
        success &= check_link_error(id, on_link_error);
    }

    glDeleteShader(compute);
    gl::gl_check();

    if (!success)
    {
        id = 0;
    }

    return success;
}

auto gl::basic_shader::validate_uniform_set(const std::string& name) -> int
{
    expect(get_id() != 0, "current shader must be valid");
    int loc = glGetUniformLocation(get_id(), name.c_str());
    gl::gl_check();
    expect(loc != -1, fmt::format("invalid name for uniform {}", name));
    return loc;
}

auto gl::basic_shader::set_int(const std::string& name, int value) -> gl::basic_shader&
{
    glProgramUniform1i(get_id(), validate_uniform_set(name), value);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_int(const std::string& name, glm::vec<2, std::int32_t> value) -> gl::basic_shader&
{
    glProgramUniform2i(get_id(), validate_uniform_set(name), value[0], value[1]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_int(const std::string& name, glm::vec<3, std::int32_t> value) -> gl::basic_shader&
{
    glProgramUniform3i(get_id(), validate_uniform_set(name), value[0], value[1], value[2]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_int(const std::string& name, glm::vec<4, std::int32_t> value) -> gl::basic_shader&
{
    glProgramUniform4i(get_id(), validate_uniform_set(name), value[0], value[1], value[2], value[3]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_float(const std::string& name, float value) -> gl::basic_shader&
{
    glProgramUniform1f(get_id(), validate_uniform_set(name), value);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_float(const std::string& name, glm::vec2 value) -> gl::basic_shader&
{
    glProgramUniform2f(get_id(), validate_uniform_set(name), value[0], value[1]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_float(const std::string& name, glm::vec3 value) -> gl::basic_shader&
{
    glProgramUniform3f(get_id(), validate_uniform_set(name), value[0], value[1], value[2]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_float(const std::string& name, glm::vec4 value) -> gl::basic_shader&
{
    glProgramUniform4f(get_id(), validate_uniform_set(name), value[0], value[1], value[2], value[3]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_uint(const std::string& name, uint value) -> gl::basic_shader&
{
    glProgramUniform1ui(get_id(), validate_uniform_set(name), value);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_uint(const std::string& name, glm::vec<2, std::uint32_t> value) -> gl::basic_shader&
{
    glProgramUniform2ui(get_id(), validate_uniform_set(name), value[0], value[1]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_uint(const std::string& name, glm::vec<3, std::uint32_t> value) -> gl::basic_shader&
{
    glProgramUniform3ui(get_id(), validate_uniform_set(name), value[0], value[1], value[2]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_vec_uint(const std::string& name, glm::vec<4, std::uint32_t> value) -> gl::basic_shader&
{
    glProgramUniform4ui(get_id(), validate_uniform_set(name), value[0], value[1], value[2], value[3]);
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_matrix(const std::string& name, const glm::mat2& matrix) -> gl::basic_shader&
{
    glProgramUniformMatrix2fv(get_id(), validate_uniform_set(name), 1, GL_FALSE, glm::value_ptr(matrix));
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_matrix(const std::string& name, const glm::mat3& matrix) -> gl::basic_shader&
{
    glProgramUniformMatrix3fv(get_id(), validate_uniform_set(name), 1, GL_FALSE, glm::value_ptr(matrix));
    gl::gl_check();
    return *this;
}

auto gl::basic_shader::set_matrix(const std::string& name, const glm::mat4& matrix) -> gl::basic_shader&
{
    glProgramUniformMatrix4fv(get_id(), validate_uniform_set(name), 1, GL_FALSE, glm::value_ptr(matrix));
    gl::gl_check();
    return *this;
}

