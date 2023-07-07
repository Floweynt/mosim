#pragma once
#include <GL/glew.h>
//
#include <GL/gl.h>
#include <concepts>
#include <functional>
#include <iostream>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#define _gl_resource_object(T)                                                                                                                       \
    T(const T&) = delete;                                                                                                                            \
    T& operator=(const T&) = delete;                                                                                                                 \
    T(T&& rhs) noexcept                                                                                                                                \
    {                                                                                                                                                \
        id = rhs.id;                                                                                                                                   \
        rhs.id = 0;                                                                                                                                    \
    }                                                                                                                                                \
    T& operator=(T&& rhs) noexcept                                                                                                                     \
    {                                                                                                                                                \
        id = rhs.id;                                                                                                                                   \
        rhs.id = 0;                                                                                                                                    \
        return *this;                                                                                                                                \
    }

namespace gl
{
    class gl_error : public std::runtime_error
    {
        int id;

    public:
        gl_error(int err_id, const char* msg) : std::runtime_error(msg), id(err_id) {}
        [[nodiscard]] constexpr auto get_id() const { return id; }
    };

    inline void gl_check()
    {
#ifdef DEBUG
        auto error = glGetError();
        if (error == 0)
            return;

        static constexpr const char* ERROR_MSG[] = {
            "GL_INVALID_ENUM(1280)",
            "GL_INVALID_VALUE(1281)",
            "GL_INVALID_OPERATION(1282)",
            "GL_STACK_OVERFLOW(1283)",
            "GL_STACK_UNDERFLOW(1284)",
            "GL_OUT_OF_MEMORY(1285)",
            "GL_INVALID_FRAMEBUFFER_OPERATION(1286)",
        };

        const char* name = "unknown gl error";
        if (error >= 1280 && error <= (1280 + sizeof(ERROR_MSG) / sizeof(ERROR_MSG[0])))
            name = ERROR_MSG[error - 1280];

        // std::cerr << ERROR_MSG[error - 1280] << '\n';
        throw gl_error(error, name);
#endif
    }

    class buffer_object
    {
        GLuint id{};

    public:
        buffer_object()
        {
            glCreateBuffers(1, &id);
            gl_check();
        }

        ~buffer_object()
        {
            if (id)
            {
                glDeleteBuffers(1, &id);
                gl_check();
            }
        }

        _gl_resource_object(buffer_object);
        // initialize

        void init(size_t size, void* ptr, GLenum usage) const
        {
            glNamedBufferData(id, size, ptr, usage);
            gl_check();
        }

        template <typename T>
        void init(const std::span<T>& data, GLenum usage)
        {
            init(data.size_bytes(), (void*)data.data(), usage);
        }

        template <typename T>
            requires(std::is_standard_layout_v<T>)
        void init(const T& object, GLenum usage)
        {
            init(sizeof(T), &object, usage);
        }

        void init_empty(size_t size, GLenum usage) const { init(size, nullptr, usage); }

        // write operations
        void write(size_t offset, size_t size, const void* ptr) const
        {
            glNamedBufferSubData(id, offset, size, ptr);
            gl_check();
        }

        template <typename T>
        void write(size_t offset, const std::span<T>& data)
        {
            write(offset, data.size_bytes(), (void*)data.data());
        }

        template <typename T>
            requires(std::is_standard_layout_v<T>)
        void write(size_t offset, const T& object)
        {
            write(offset, sizeof(T), (void*)&object);
        }

        // read operations
        void read_to(size_t offset, size_t size, void* ptr) const
        {
            glGetNamedBufferSubData(id, offset, size, ptr);
            gl_check();
        }

        template <typename T>
            requires(!std::is_const_v<T>)
        void read_to(size_t offset, const std::span<T>& data)
        {
            read_to(offset, data.size_bytes(), (void*)data.data());
        }

        template <typename T>
            requires(std::is_standard_layout_v<T>)
        void read_to(size_t offset, T& object)
        {
            read_to(offset, sizeof(T), (void*)&object);
        }

        template <typename T>
            requires(std::is_arithmetic_v<T>)
        auto read(size_t offset) -> T
        {
            T inst{};
            read_to(offset, inst);
            return inst;
        }

        [[nodiscard]] constexpr auto get_id() const { return id; }
    };

    class vertex_array_object
    {
        GLuint id{};

    public:
        vertex_array_object()
        {
            glCreateVertexArrays(1, &id);
            gl_check();
        }
        ~vertex_array_object()
        {
            if (id)
            {
                glDeleteBuffers(1, &id);
                gl_check();
            }
        }

        _gl_resource_object(vertex_array_object);

        void format(uint32_t attribute, size_t elements, GLenum type, bool normalized, size_t offset) const
        {
            glVertexArrayAttribFormat(id, attribute, elements, type, normalized ? GL_TRUE : GL_FALSE, offset);
            gl_check();
        }

        void int_format(uint32_t attribute, size_t elements, GLenum type, size_t offset) const
        {
            glVertexArrayAttribIFormat(id, attribute, elements, type, offset);
            gl_check();
        }

        void long_format(uint32_t attribute, size_t elements, GLenum type, size_t offset) const
        {
            glVertexArrayAttribLFormat(id, attribute, elements, type, offset);
            gl_check();
        }

        void bind_vbo(const buffer_object& object, uint32_t bind_point, intptr_t offset, size_t stride) const
        {
            glVertexArrayVertexBuffer(id, bind_point, object.get_id(), offset, stride);
            gl_check();
        }

        void enable_attribute(uint32_t attribute) const
        {
            glEnableVertexArrayAttrib(id, attribute);
            gl_check();
        }

        void set_attribute_buffer_bind(uint32_t attribute, uint32_t target_bind_point) const
        {
            glVertexArrayAttribBinding(id, attribute, target_bind_point);
            gl_check();
        }

        void disable_attribute(uint32_t attribute) const
        {
            glDisableVertexArrayAttrib(id, attribute);
            gl_check();
        }

        void draw(GLenum primitive, size_t first, size_t count) const
        {
            glBindVertexArray(id);
            gl_check();
            glDrawArrays(primitive, first, count);
            gl_check();
            glBindVertexArray(0);
            gl_check();
        }

        [[nodiscard]] constexpr auto get_id() const { return id; }
    };

    class texture_unit
    {
        uint32_t id;

    public:
        constexpr texture_unit(uint32_t id) : id(id){};
        [[nodiscard]] constexpr auto get_unit() const { return GL_TEXTURE0 + id; };
        [[nodiscard]] constexpr auto get_id() const { return id; }
    };

    inline static constexpr texture_unit TEXTURE_UNIT_0 = texture_unit(0);
    inline static constexpr texture_unit TEXTURE_UNIT_1 = texture_unit(1);
    inline static constexpr texture_unit TEXTURE_UNIT_2 = texture_unit(2);
    inline static constexpr texture_unit TEXTURE_UNIT_3 = texture_unit(3);
    inline static constexpr texture_unit TEXTURE_UNIT_4 = texture_unit(4);
    inline static constexpr texture_unit TEXTURE_UNIT_5 = texture_unit(5);
    inline static constexpr texture_unit TEXTURE_UNIT_6 = texture_unit(6);
    inline static constexpr texture_unit TEXTURE_UNIT_7 = texture_unit(7);

    class texture_object
    {
    private:
        GLuint id{};
        const texture_unit& default_unit;

    public:
        texture_object(GLenum type, const texture_unit& default_unit = TEXTURE_UNIT_0) : default_unit(default_unit)
        {
            glCreateTextures(type, 1, &id);
            gl_check();
        }
        ~texture_object()
        {
            if (id)
            {
                glDeleteTextures(1, &id);
                gl_check();
            }
        }

        void bind() const { bind(default_unit); }
        void bind(const texture_unit& unit) const
        {
            glBindTextureUnit(unit.get_id(), id);
            gl_check();
        }

        void set_float(GLenum pname, GLfloat param) const
        {
            glTextureParameterf(id, pname, param);
            gl_check();
        }

        void set_int(GLenum pname, GLint param) const
        {
            glTextureParameteri(id, pname, param);
            gl_check();
        }

        void set_storage_format(size_t levels, GLenum internalformat, size_t width, size_t height) const
        {
            glTextureStorage2D(id, levels, internalformat, width, height);
            gl_check();
        }

        void upload_texture(int32_t level, int32_t xoffset, int32_t yoffset, size_t width, size_t height, GLenum format, GLenum type,
                            const void* pixels) const
        {
            glTextureSubImage2D(id, level, xoffset, yoffset, width, height, format, type, pixels);
            gl_check();
        }

        [[nodiscard]] constexpr auto get_id() const { return id; }
        [[nodiscard]] constexpr auto get_default_unit() const -> const auto& { return default_unit; }
    };
} // namespace gl
