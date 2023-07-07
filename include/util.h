#pragma once

#include <stdexcept>
#include <string>

constexpr void expect_true(bool cond, const std::string& str)
{
    if (cond != true)
        throw std::runtime_error(str);
}

constexpr void expect_false(bool cond, const std::string& str)
{
    if (cond != false)
        throw std::runtime_error(str);
}

constexpr void expect(bool cond, const std::string& str) { expect_true(cond, str); }

template <typename T>
constexpr auto decay_val(const T& val) -> std::decay_t<T>
{
    return std::decay_t<T>(val);
}
