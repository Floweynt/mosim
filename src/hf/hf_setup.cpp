// cSpell:ignore fgamma dvec
#include "hf/basis.h"
#include "hf/chem.h"
#include "hf/gamma.h"
#include "hf/math_util.h"
#include <cstdint>
#include <glm/ext.hpp>
#include <glm/gtx/norm.hpp>
#include <vector>

#define FGAMMA_CACHE_FOR(_size, expr)                                                                                                                \
    std::vector<double> _fgamma_lookup_table((_size) + 1);                                                                                           \
    for (size_t i = 0; i < _fgamma_lookup_table.size(); i++)                                                                                         \
    {                                                                                                                                                \
        _fgamma_lookup_table[i] = Fgamma(i, expr);                                                                                                   \
    } 

#define FGAMMA_CACHE_GET(index) _fgamma_lookup_table[index]

namespace hf
{
    INLINE auto b0(int i, int r, double g) -> double { return fact_ratio2(i, r) * pow_i(4 * g, r - i); }

    INLINE auto force_b(int64_t i, int64_t l1, int64_t l2, double p, double a, double b, int64_t r, double g) -> double
    {
        return binomial_prefactor(i, l1, l2, p - a, p - b) * b0(i, r, g);
    }

    INLINE auto b_term(int64_t i1, int64_t i2, int64_t r1, int64_t r2, int64_t u, int64_t l1, int64_t l2, int64_t l3, int64_t l4, double px,
        double ax, double bx, double qx, double cx, double dx, double gamma1, double gamma2, double delta) -> double
    {
        return force_b(i1, l1, l2, px, ax, bx, r1, gamma1) * pow_i(-1, i2) * force_b(i2, l3, l4, qx, cx, dx, r2, gamma2) * pow_i(-1, u) *
               fact_ratio2(i1 + i2 - 2 * (r1 + r2), u) * pow_i(qx - px, i1 + i2 - 2 * (r1 + r2) - 2 * u) / pow_i(delta, i1 + i2 - 2 * (r1 + r2) - u);
    }

    INLINE auto b_array(int64_t l1, int64_t l2, int64_t l3, int64_t l4, double p, double a, double b, double q, double c, double d, double g1,
        double g2, double delta) -> std::vector<double>
    {
        int max_idx = l1 + l2 + l3 + l4 + 1;
        std::vector<double> arr(max_idx, 0);
        for (int i1 = 0; i1 < l1 + l2 + 1; i1++)
        {
            for (int i2 = 0; i2 < l3 + l4 + 1; i2++)
            {
                for (int r1 = 0; r1 < i1 / 2 + 1; r1++)
                {
                    for (int r2 = 0; r2 < i2 / 2 + 1; r2++)
                    {
                        for (int u = 0; u < (i1 + i2) / 2 - r1 - r2 + 1; u++)
                        {
                            int i = i1 + i2 - 2 * (r1 + r2) - u;
                            arr[i] += b_term(i1, i2, r1, r2, u, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta);
                        }
                    }
                }
            }
        }
        return arr;
    }

    INLINE auto a_term(int64_t i, int64_t r, int64_t u, int64_t l1, int64_t l2, double pax, double pbx, double cpx, double gamma) -> double
    {
        return pow_i(-1, i) * binomial_prefactor(i, l1, l2, pax, pbx) * pow_i(-1, u) * fact(i) * pow_i(cpx, i - 2 * r - 2 * u) *
               pow_i(0.25 / gamma, r + u) / fact(r) / fact(u) / fact(i - 2 * r - 2 * u);
    }

    INLINE auto a_array(int64_t l1, int64_t l2, double pa, double pb, double cp, double g) -> std::vector<double>
    {
        int max_idx = l1 + l2 + 1;
        std::vector<double> arr(max_idx, 0);
        for (int i = 0; i < max_idx; i++)
        {
            for (int r = 0; r <= i / 2; r++)
            {
                for (int u = 0; u <= (i - 2 * r) / 2; u++)
                {
                    int iI = i - 2 * r - u;
                    arr[iI] += a_term(i, r, u, l1, l2, pa, pb, cp, g);
                }
            }
        }
        return arr;
    }

    INLINE auto overlap_1d(int64_t l1, int64_t l2, double x1, double x2, double gamma) -> double
    {
        double sum = 0;
        for (int i = 0; i < (1 + floor(0.5 * (l1 + l2))); i++)
        {
            sum += binomial_prefactor(2 * i, l1, l2, x1, x2) * fact2(2 * i - 1) / pow_i(2 * gamma, i);
        }
        return sum;
    }

    INLINE auto overlap(double alpha1, int l1, int m1, int n1, const glm::dvec3& a, double alpha2, int l2, int m2, int n2, const glm::dvec3& b)
        -> double
    {
        double rab2 = glm::distance2(a, b);
        double gamma = alpha1 + alpha2;
        glm::dvec3 p = gaussian_product(alpha1, a, alpha2, b);
        double pre = pow(M_PI / gamma, 1.5) * exp(-alpha1 * alpha2 * rab2 / gamma);
        double wx = overlap_1d(l1, l2, p.x - a.x, p.x - b.x, gamma);
        double wy = overlap_1d(m1, m2, p.y - a.y, p.y - b.y, gamma);
        double wz = overlap_1d(n1, n2, p.z - a.z, p.z - b.z, gamma);
        return pre * wx * wy * wz;
    }

    INLINE auto gto_overlap(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
    {
        return overlap(
            gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(), pos2);
    }

    auto cgf_overlap(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2) -> double
    {
        double sum = 0;

        for (const auto& gto1 : cgf1.orbs())
        {
            for (const auto& gto2 : cgf2.orbs())
            {
                sum += gto1.factor() * gto2.factor() * gto_overlap(gto1, gto2, cgf1.get_pos(), cgf2.get_pos());
            }
        }

        return sum;
    }

    INLINE auto gto_kinetic(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
    {
        double term0 = gto2.get_alpha() * (2 * (gto2.get_l() + gto2.get_m() + gto2.get_n()) + 3) * gto_overlap(gto1, gto2, pos1, pos2);
        double term1 = -2 * pow_i(gto2.get_alpha(), 2) *
                       (overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l() + 2, gto2.get_m(),
                            gto2.get_n(), pos2) +
                           overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(), gto2.get_m() + 2,
                               gto2.get_n(), pos2) +
                           overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(), gto2.get_m(),
                               gto2.get_n() + 2, pos2));
        double term2 = -0.5 * (gto2.get_l() * (gto2.get_l() - 1) *
                                      overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l() - 2,
                                          gto2.get_m(), gto2.get_n(), pos2) +
                                  gto2.get_m() * (gto2.get_m() - 1) *
                                      overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(),
                                          gto2.get_m() - 2, gto2.get_n(), pos2) +
                                  gto2.get_n() * (gto2.get_n() - 1) *
                                      overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(),
                                          gto2.get_m(), gto2.get_n() - 2, pos2));
        return term0 + term1 + term2;
    }

    auto cgf_kinetic(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2) -> double
    {
        double sum = 0;
        for (const auto& gto1 : cgf1.orbs())
        {
            for (const auto& gto2 : cgf2.orbs())
            {
                sum += gto1.factor() * gto2.factor() * gto_kinetic(gto1, gto2, cgf1.get_pos(), cgf2.get_pos());
            }
        }

        return sum;
    }

    INLINE auto nuclear(const glm::dvec3& a, double norm1, int l1, int m1, int n1, double alpha1, const glm::dvec3& b, double norm2, int l2, int m2,
        int n2, double alpha2, const glm::dvec3& c) -> double
    {
        double gamma = alpha1 + alpha2;
        glm::dvec3 p = gaussian_product(alpha1, a, alpha2, b);
        double rab2 = glm::distance2(a, b);
        double rcp2 = glm::distance2(c, p);
        std::vector<double> ax = a_array(l1, l2, p.x - a.x, p.x - b.x, p.x - c.x, gamma);
        std::vector<double> ay = a_array(m1, m2, p.y - a.y, p.y - b.y, p.y - c.y, gamma);
        std::vector<double> az = a_array(n1, n2, p.z - a.z, p.z - b.z, p.z - c.z, gamma);
        double sum = 0.0;

        FGAMMA_CACHE_FOR(l1 + l2 + m1 + m2 + n1 + n2, rcp2 * gamma);

        for (int i = 0; i <= l1 + l2; i++)
        {
            for (int j = 0; j <= m1 + m2; j++)
            {
                for (int k = 0; k <= n1 + n2; k++)
                {
                    sum += ax[i] * ay[j] * az[k] * FGAMMA_CACHE_GET(i + j + k);
                }
            }
        }

        return -norm1 * norm2 * 2 * M_PI / gamma * exp(-alpha1 * alpha2 * rab2 / gamma) * sum;
    }

    INLINE auto gto_nuclear(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& c, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
    {
        return nuclear(pos1, gto1.get_norm(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(), 
                      pos2, gto2.get_norm(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), 
                      c);
    }

    auto cgf_nuclear(contracted_gaussian_functions& cgf1, contracted_gaussian_functions& cgf2, const atom& a) -> double
    {
        double sum = 0;
#pragma omp parallel for collapse(2) reduction(+ : sum)
        for (const auto& gto1 : cgf1.orbs())
        {
            for (const auto& gto2 : cgf2.orbs())
            {
                sum += gto1.get_coeff() * gto2.get_coeff() * gto_nuclear(gto1, gto2, a.pos(), cgf1.get_pos(), cgf2.get_pos());
            }
        }

        return sum * a.proton_count();
    }

    INLINE auto coulomb_repulsion(const glm::dvec3& a, double norma, int64_t la, int64_t ma, int64_t na, double alphaa, const glm::dvec3& b,
        double normb, int64_t lb, int64_t mb, int64_t nb, double alphab, const glm::dvec3& c, double normc, int64_t lc, int64_t mc, int64_t nc,
        double alphac, const glm::dvec3& d, double normd, int64_t ld, int64_t md, int64_t nd, double alphad) -> double
    {
        double rab2 = glm::distance2(a, b);
        double rcd2 = glm::distance2(c, d);
        glm::dvec3 p = gaussian_product(alphaa, a, alphab, b);
        glm::dvec3 q = gaussian_product(alphac, c, alphad, d);
        double rpq2 = glm::distance2(p, q);
        double gamma1 = alphaa + alphab;
        double gamma2 = alphac + alphad;
        double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);
        std::vector<double> bx = b_array(la, lb, lc, ld, p.x, a.x, b.x, q.x, c.x, d.x, gamma1, gamma2, delta);
        std::vector<double> by = b_array(ma, mb, mc, md, p.y, a.y, b.y, q.y, c.y, d.y, gamma1, gamma2, delta);
        std::vector<double> bz = b_array(na, nb, nc, nd, p.z, a.z, b.z, q.z, c.z, d.z, gamma1, gamma2, delta);
        double sum = 0.0;

        FGAMMA_CACHE_FOR(la + lb + lc + ld + ma + mb + mc + md + na + nb + nc + nd, 0.25 * rpq2 / delta);

        for (int i = 0; i <= (la + lb + lc + ld); i++)
        {
            for (int j = 0; j <= (ma + mb + mc + md); j++)
            {
                for (int k = 0; k <= (na + nb + nc + nd); k++)
                {
                    sum += bx[i] * by[j] * bz[k] * FGAMMA_CACHE_GET(i + j + k);
                }
            }
        }

        static constexpr auto SQRT_PI_FIFTH = 17.493418327624862;

        return 2 * SQRT_PI_FIFTH / (gamma1 * gamma2 * sqrt(gamma1 + gamma2)) * exp(-alphaa * alphab * rab2 / gamma1) *
               exp(-alphac * alphad * rcd2 / gamma2) * sum * norma * normb * normc * normd;
    }

    INLINE auto gto_repulsion(const gto_data& gto1, const gto_data& gto2, const gto_data& gto3, const gto_data& gto4, const glm::dvec3& pos1,
        const glm::dvec3& pos2, const glm::dvec3& pos3, const glm::dvec3& pos4) -> double
    {
        return coulomb_repulsion(pos1, gto1.get_norm(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(), pos2, gto2.get_norm(),
            gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), pos3, gto3.get_norm(), gto3.get_l(), gto3.get_m(), gto3.get_n(),
            gto3.get_alpha(), pos4, gto4.get_norm(), gto4.get_l(), gto4.get_m(), gto4.get_n(), gto4.get_alpha());
    }

    auto cgf_repulsion(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2,
        const contracted_gaussian_functions& cgf3, const contracted_gaussian_functions& cgf4) -> double
    {
        double sum = 0;
        int size1 = cgf1.orbs().size();
        int size2 = cgf2.orbs().size();
        int size3 = cgf3.orbs().size();
        int size4 = cgf4.orbs().size();
#pragma omp parallel for collapse(4) reduction(+ : sum)
        for (int i = 0; i < size1; i++)
        {
            for (int j = 0; j < size2; j++)
            {
                for (int k = 0; k < size3; k++)
                {
                    for (int l = 0; l < size4; l++)
                    {
                        double pre =
                            cgf1.orbs()[i].get_coeff() * cgf2.orbs()[j].get_coeff() * cgf3.orbs()[k].get_coeff() * cgf4.orbs()[l].get_coeff();
                        sum += pre * gto_repulsion(cgf1.orbs()[i], cgf2.orbs()[j], cgf3.orbs()[k], cgf4.orbs()[l], cgf1.get_pos(), cgf2.get_pos(),
                                         cgf3.get_pos(), cgf4.get_pos());
                    }
                }
            }
        }
        return sum;
    }
} // namespace hf
