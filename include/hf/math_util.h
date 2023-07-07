#pragma once
#include <glm/glm.hpp>

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

    int res = 1;
    while (val >= 11)
    {
        res *= val--;
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

    int res = 1;
    while (val >= 11)
    {
        res *= val;
        val -= 2;
    }

    return res * TABLE[val];
}

constexpr auto binomial(int a, int b) -> int { return fact(a) / (fact(b) * fact(a - b)); }

constexpr auto gaussian_product_center(const double& alpha1, const glm::dvec3& a, const double& alpha2, const glm::dvec3& b) -> glm::dvec3
{
    double gamma = alpha1 + alpha2;
    return {(alpha1 * a.x + alpha2 * b.x) / gamma, (alpha1 * a.y + alpha2 * b.y) / gamma, (alpha1 * a.z + alpha2 * b.z) / gamma};
}

constexpr auto binomial_prefactor(int s, int ia, int ib, double xpa, double xpb) -> double
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
