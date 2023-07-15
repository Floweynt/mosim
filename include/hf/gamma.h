// cSpell:disable
#pragma once

#include <cmath>
#include <iostream>
#include <vector>

namespace detail
{
    constexpr auto gcf(double a, double x) -> double
    {
        constexpr double FP_MIN = std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();
        double gln = std::lgamma(a);

        double b = x + 1.0 - a;
        double c = 1.0 / FP_MIN;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1;; i++)
        {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::fabs(d) < FP_MIN)
            {
                d = FP_MIN;
            }
            c = b + an / c;
            if (std::fabs(c) < FP_MIN)
            {
                c = FP_MIN;
            }
            d = 1.0 / d;
            double del = d * c;
            h *= del;
                                                                                   // :P
            if (std::fabs(del - 1.0) <= /*std::numeric_limits<double>::epsilon()*/ 0.00001)
            {
                break;
            }
        }
        return std::exp(-x + a * std::log(x) - gln) * h;
    }

    constexpr auto gammp_approx(double a, double x, int psig) -> double
    {
        constexpr int NGAU = 18;
        constexpr double y[] = {0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421, 0.082502225484340941,
                                0.12007019910960293,   0.16415283300752470,  0.21442376986779355,  0.27051082840644336,  0.33199876341447887,
                                0.39843234186401943,   0.46931971407375483,  0.54413605556657973,  0.62232745288031077,  0.70331500465597174,
                                0.78649910768313447,   0.87126389619061517,  0.95698180152629142};
        constexpr double w[] = {0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734, 0.034213810770299537,
                                0.040875750923643261,  0.047235083490265582, 0.053244713977759692, 0.058860144245324798, 0.064039797355015485,
                                0.068745323835736408,  0.072941885005653087, 0.076598410645870640, 0.079687828912071670, 0.082187266704339706,
                                0.084078218979661945,  0.085346685739338721, 0.085983275670394821};
        double xu;
        double a1 = a - 1.0;
        double lna1 = log(a1);
        double sqrta1 = sqrt(a1);
        double gln = std::lgamma(a);
        if (x > a1)
        {
            xu = std::max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
        }
        else
        {
            xu = std::max(0., std::min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
        }

        double sum = 0;

        for (int j = 0; j < NGAU; j++)
        {
            double t = x + (xu - x) * y[j];
            sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
        }

        double ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
        return (psig ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans));
    }

    constexpr auto gser(double a, double x) -> double
    {
        double sum = 0;
        double gln = std::lgamma(a);
        double ap = a;
        double del = sum = 1.0 / a;
        while (true)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * std::numeric_limits<double>::epsilon())
            {
                return sum * std::exp(-x + a * std::log(x) - gln);
            }
        }
    }

    constexpr auto gammp(double a, double x) -> double
    {
        constexpr int ASWITCH = 100;
        if (x < 0.0 || a <= 0.0)
        {
            throw std::domain_error("Bad value in Fgamma!");
            return 0.0;
        }
        if (x == 0.0)
        {
            return 0.0;
        }
        if (a >= ASWITCH)
        {
            return gammp_approx(a, x, 1);
        }
        if (x < a + 1.0)
        {
            return gser(a, x);
        }

        return 1.0 - gcf(a, x);
    }

    constexpr auto gamm_inc(double a, double x) -> double
    {
        double gammap = gammp(a, x);
        return std::tgamma(a) * gammap;
    }
} // namespace detail

constexpr auto Fgamma(int64_t m, double x) -> double
{
    constexpr double tiny = 0.00000001;
    x = std::max(std::abs(x), tiny);
    return 0.5 * pow(x, -m - 0.5) * detail::gamm_inc(m + 0.5, x);
}
