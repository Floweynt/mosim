#pragma once
#include <array>
#include <build_config.h>
#include <glm/glm.hpp>
#include <type_traits>

INLINE
#if __cplusplus > 202002L
constexpr
#endif
    auto
    fact(int64_t val) -> int64_t
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

#if __cplusplus > 202002L
    if consteval
    {
        int64_t res = 1;
        while (val >= 11)
        {
            res *= val--;
        }

        return res * TABLE[val];
    }
#endif

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

INLINE
#if __cplusplus > 202002L
constexpr
#endif
    auto
    fact2(int64_t val) -> int64_t
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

#if __cplusplus > 202002L
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
#endif

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

INLINE
#if __cplusplus > 202002L
constexpr
#endif

    auto
    binomial(int64_t a, int64_t b) -> int64_t
{
#if __cplusplus > 202002L
    if consteval
    {
        return fact(a) / (fact(b) * fact(a - b));
    }
#endif

    static constexpr size_t CACHE_SIZE = 128;
    static std::array<int64_t, CACHE_SIZE * CACHE_SIZE> cache{};

    bool cachable = (0 <= a && a < CACHE_SIZE) && (0 <= b && b < CACHE_SIZE);
    if (cachable && cache[a * CACHE_SIZE + b])
    {
        return cache[a * CACHE_SIZE + b];
    }

    if (cachable)
    {
        return cache[a * CACHE_SIZE + b] = fact(a) / (fact(b) * fact(a - b));
    }
    return fact(a) / (fact(b) * fact(a - b));
}

INLINE
#if __cplusplus > 202002L
constexpr
#endif

    auto
    gaussian_product(const double& alpha1, const glm::dvec3& a, const double& alpha2, const glm::dvec3& b) -> glm::dvec3
{
    double gamma = alpha1 + alpha2;
    return {(alpha1 * a.x + alpha2 * b.x) / gamma, (alpha1 * a.y + alpha2 * b.y) / gamma, (alpha1 * a.z + alpha2 * b.z) / gamma};
}

template <typename T>
    requires(std::is_integral_v<T>)
INLINE
#if __cplusplus > 202002L
    constexpr
#endif
    auto
    pow_i(double base, T exp)
{
    double result = 1;

    if (exp < 0)
    {
        base = 1 / base;
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

INLINE
#if __cplusplus > 202002L
constexpr
#endif

    auto
    binomial_prefactor(int64_t s, int64_t ia, int64_t ib, double xpa, double xpb) -> double
{
    double sum = 0.;
    for (int64_t t = 0; t < s + 1; t++)
    {
        if ((s - ia <= t) && (t <= ib))
        {
            sum += binomial(ia, s - t) * binomial(ib, t) * pow_i(xpa, ia - s + t) * pow_i(xpb, ib - t);
        }
    }
    return sum;
}

#if __cplusplus > 202002L
constexpr
#endif
    INLINE auto
    fact_ratio2(int64_t a, int64_t b) -> double
{
#if __cplusplus > 202002L
    if consteval
    {
        return fact(a) / fact(b) / fact(a - 2 * b);
    }
#endif

    static constexpr size_t CACHE_SIZE = 128;
    static std::array<int64_t, CACHE_SIZE * CACHE_SIZE> cache{};

    bool cachable = (0 <= a && a < CACHE_SIZE) && (0 <= b && b < CACHE_SIZE);
    if (cachable && cache[a * CACHE_SIZE + b])
    {
        return cache[a * CACHE_SIZE + b];
    }

    if (cachable)
    {
        return cache[a * CACHE_SIZE + b] = fact(a) / fact(b) / fact(a - 2 * b);
    }
    return fact(a) / fact(b) / fact(a - 2 * b);
}
