#pragma once
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <vector>

namespace gl
{
    using namespace glm;
    class matrix_stack
    {
        std::vector<mat4> st;
        mat4 current;

    public:
        constexpr matrix_stack() : current(glm::identity<mat4>()) {}

        // pushes the current matrix, then multiples by mat
         void push(const mat4& mat)
        {
            st.push_back(current);
            current = current * mat;
        }

        constexpr auto push() -> auto&
        {
            st.push_back(current);
            return *this;
        }

        constexpr void pop()
        {
            current = st.back();
            st.pop_back();
        }

         void mult(const mat4& mat) { current *= mat; }

        [[nodiscard]] constexpr auto curr() const -> const auto& { return current; }

         auto rotate(float theta, const vec3& axis) -> auto&
        {
            current *= glm::rotate(current, theta, axis);
            return *this;
        }

         auto rotate(const glm::quat& quat) -> auto&
        {
            current *= glm::toMat4(quat);
            return *this;
        }

         auto translate(const vec3& offset) -> auto&
        {
            current *= glm::translate(current, offset);
            return *this;
        }

         auto scale(const vec3& val) -> auto&
        {
            current *= glm::scale(current, val);
            return *this;
        }

         auto rotate(float theta, float x, float y, float z) -> auto& { return rotate(theta, {x, y, z}); }
         auto translate(float x, float y, float z) -> auto& { return translate({x, y, z}); }
         auto scale(float x, float y, float z) -> auto& { return scale({x, y, z}); }
         auto scale(float amount) -> auto& { return scale(amount, amount, amount); }
    };
} // namespace gl

