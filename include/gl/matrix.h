#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/quaternion_transform.hpp>
#include <glm/fwd.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <vector>

namespace gl
{
    using namespace glm;
    class matrix_stack
    {
        std::vector<mat4> st;
        mat4 current;

    public:
        constexpr matrix_stack() : st(), current(glm::identity<mat4>()) {}

        // pushes the current matrix, then multiples by mat
        constexpr void push(const mat4& mat)
        {
            st.push_back(current);
            current = current * mat;
        }

        constexpr auto& push()
        {
            st.push_back(current);
            return *this;
        }

        constexpr void pop()
        {
            current = st.back();
            st.pop_back();
        }

        constexpr void mult(const mat4& mat) { current *= mat; }

        [[nodiscard]] constexpr auto curr() const -> const auto& { return current; }

        constexpr auto rotate(float theta, const vec3& axis) -> auto&
        {
            current *= glm::rotate(current, theta, axis);
            return *this;
        }

        constexpr auto rotate(const glm::quat& quat) -> auto&
        {
            current *= glm::toMat4(quat);
            return *this;
        }

        constexpr auto translate(const vec3& offset) -> auto&
        {
            current *= glm::translate(current, offset);
            return *this;
        }

        constexpr auto scale(const vec3& val) -> auto&
        {
            current *= glm::scale(current, val);
            return *this;
        }

        constexpr auto rotate(float theta, float x, float y, float z) -> auto& { return rotate(theta, {x, y, z}); }
        constexpr auto translate(float x, float y, float z) -> auto& { return translate({x, y, z}); }
        constexpr auto scale(float x, float y, float z) -> auto& { return scale({x, y, z}); }
        constexpr auto scale(float amount) -> auto& { return scale(amount, amount, amount); }
    };
} // namespace gl

#endif
