#pragma once

#include <cstdint>
#include <functional>
#include <glm/ext.hpp>
#include <optional>
#include <stdexcept>
#include <string>

namespace gl
{
    class shader_error : public std::runtime_error
    {
        std::string log;
        inline shader_error(const std::string& msg, std::string log) : runtime_error(msg), log(std::move(log)) {}
        friend class basic_shader;

    public:
        [[nodiscard]] constexpr auto get_log() const -> const auto& { return log; }
    };

    using i32vec2 = glm::vec<2, std::int32_t>;
    using i32vec3 = glm::vec<3, std::int32_t>;
    using i32vec4 = glm::vec<4, std::int32_t>;
    using u32vec2 = glm::vec<2, std::uint32_t>;
    using u32vec3 = glm::vec<3, std::uint32_t>;
    using u32vec4 = glm::vec<4, std::uint32_t>;

    class basic_shader
    {
        auto validate_uniform_set(const std::string& name) -> int;
    protected:
        unsigned int id{0};
        using err_hdl_t = const std::function<void(std::string)>&;

        virtual auto do_compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool = 0;
    public:
        constexpr basic_shader() = default;

        auto compile() -> bool;
        auto compile(err_hdl_t on_error) -> bool;
        auto compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool;
        void use() const;
        [[nodiscard]] constexpr auto get_id() const { return id; }
        void dispose();
        
        auto set_int(const std::string& name, std::int32_t value) -> basic_shader&;
        auto set_vec_int(const std::string& name, i32vec2 value) -> basic_shader&;
        auto set_vec_int(const std::string& name, i32vec3 value) -> basic_shader&;
        auto set_vec_int(const std::string& name, i32vec4 value) -> basic_shader&;

        auto set_float(const std::string& name, float value) -> basic_shader&;
        auto set_vec_float(const std::string& name, glm::vec2 value) -> basic_shader&;
        auto set_vec_float(const std::string& name, glm::vec3 value) -> basic_shader&;
        auto set_vec_float(const std::string& name, glm::vec4 value) -> basic_shader&;

        auto set_uint(const std::string& name, std::uint32_t value) -> basic_shader&;
        auto set_vec_uint(const std::string& name, u32vec2 value) -> basic_shader&;
        auto set_vec_uint(const std::string& name, u32vec3 value) -> basic_shader&;
        auto set_vec_uint(const std::string& name, u32vec4 value) -> basic_shader&;

        auto set_matrix(const std::string& name, const glm::mat2& matrix) -> basic_shader&;
        auto set_matrix(const std::string& name, const glm::mat3& matrix) -> basic_shader&;
        auto set_matrix(const std::string& name, const glm::mat4& matrix) -> basic_shader&;

        virtual ~basic_shader();
    };

    class shader : public basic_shader
    {
        std::string vertex_path;
        std::string fragment_path;
        std::string geom_path;
    protected:

        auto do_compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool override;
    public:
        constexpr shader(std::string vertex_path, std::string fragment_path, std::string geom_path = "")
            : vertex_path(std::move(vertex_path)), fragment_path(std::move(fragment_path)), geom_path(std::move(geom_path))
        {
        }

        ~shader() override = default;
    };

    class compute_shader  : public basic_shader
    {
        std::string path;
        protected:

        auto do_compile(err_hdl_t on_compile_error, err_hdl_t on_link_error) -> bool override;

    public:
        constexpr compute_shader(std::string path) : path(std::move(path)) {}

        ~compute_shader() override = default;
    };
} // namespace gl
