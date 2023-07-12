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

template <typename T>
constexpr auto swap_smaller_first(T& a, T& b)
{
    if (a > b)
    {
        std::swap(a, b);
    }
}

template <typename T>
constexpr auto swap_larger_first(T& a, T& b)
{
    if (a < b)
    {
        std::swap(a, b);
    }
}
