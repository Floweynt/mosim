#pragma once
#include <array>
#include <glm/glm.hpp>
#include <type_traits>

constexpr auto fact(int64_t val) -> int64_t
{
    constexpr int64_t TABLE[11] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628880};
    if (val < 2)
    {
        return 1;
    }
    if (val < 11)
    {
        return TABLE[val];
    }

    if consteval
    {
        int64_t res = 1;
        while (val >= 11)
        {
            res *= val--;
        }

        return res * TABLE[val];
    }

    static std::array<int64_t, 64> cache{};

    if (val < 64 && cache[val])
    {
        return cache[val];
    }

    int64_t res = 1;
    while (val >= 11)
    {
        res *= val--;
    }

    if (val < 64)
    {
        return cache[val] = res * TABLE[val];
    }

    return res * TABLE[val];
}

constexpr auto fact2(int64_t val) -> int64_t
{
    constexpr int64_t TABLE[11] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840};
    if (val < 2)
    {
        return 1;
    }
    if (val < 11)
    {
        return TABLE[val];
    }

    if consteval
    {
        int64_t res = 1;
        while (val >= 11)
        {
            res *= val;
            val -= 2;
        }
        return res;
    }

    static std::array<int64_t, 64> cache{};

    if (val < 64 && cache[val])
    {
        return cache[val];
    }

    int64_t res = 1;
    while (val >= 11)
    {
        res *= val;
        val -= 2;
    }

    if (val < 64)
    {
        return cache[val] = res * TABLE[val];
    }
    return res * TABLE[val];
}

constexpr auto binomial(int64_t a, int64_t b) -> int64_t
{
    if consteval
    {
        return fact(a) / (fact(b) * fact(a - b));
    }

    static constexpr size_t CACHE_SIZE = 64;
    static std::array<int64_t, CACHE_SIZE * CACHE_SIZE> cache{};
    if (0 < a && a < CACHE_SIZE && 0 < b && b < CACHE_SIZE && cache[a * CACHE_SIZE + b])
    {
        return cache[a * CACHE_SIZE + b];
    }

    if (0 < a && a < CACHE_SIZE && 0 < b && b < CACHE_SIZE)
    {
        return cache[a * CACHE_SIZE + b] = fact(a) / (fact(b) * fact(a - b));
    }
    return fact(a) / (fact(b) * fact(a - b));
}

constexpr auto gaussian_product_center(const double& alpha1, const glm::dvec3& a, const double& alpha2, const glm::dvec3& b) -> glm::dvec3
{
    double gamma = alpha1 + alpha2;
    return {(alpha1 * a.x + alpha2 * b.x) / gamma, (alpha1 * a.y + alpha2 * b.y) / gamma, (alpha1 * a.z + alpha2 * b.z) / gamma};
}

constexpr auto binomial_prefactor(int64_t s, int64_t ia, int64_t ib, double xpa, double xpb) -> double
{
    double sum = 0.;
    for (int t = 0; t < s + 1; t++)
    {
        if ((s - ia <= t) && (t <= ib))
        {
            sum += binomial(ia, s - t) * binomial(ib, t) * std::pow(xpa, ia - s + t) * std::pow(xpb, ib - t);
        }
    }
    return sum;
}

template<typename T> requires(std::is_integral_v<T>)
constexpr auto pow_i(double base, T exp)
{
    double result = 1;

    if(exp < 0)
    {
        base = 1/base;
        exp = -exp;
    }

    while (true)
    {
        if ((exp & 1) != 0)
        {
            result *= base;
        }
        exp >>= 1;
        if (exp == 0)
        {
            break;
        }
        base *= base;
    }

    return result;
}
