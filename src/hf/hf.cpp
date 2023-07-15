#include "hf/hf.h"
#include "hf/basis.h"
#include "hf/chem.h"
#include "hf/gamma.h"
#include <fmt/core.h>
#include <glm/ext.hpp>
#include <iostream>
#include <ranges>

using H_matrix = Eigen::MatrixXd;
using T_matrix = Eigen::MatrixXd;
using V_matrix = Eigen::MatrixXd;
using S_matrix = Eigen::MatrixXd;
using X_matrix = Eigen::MatrixXd;
using F_matrix = Eigen::MatrixXd;
using P_matrix = Eigen::MatrixXd;
using G_matrix = Eigen::MatrixXd;
using Cc_matrix = Eigen::MatrixXd;
using C_matrix = Eigen::MatrixXd;

namespace
{
    void compute_eigenvalues(const Eigen::MatrixXd& matrix, Eigen::MatrixXd& output, std::vector<double>& lambda)
    {
        auto size = matrix.rows();
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
        lambda.resize(size);

        for (size_t i = 0; i < size; i++)
        {
            lambda[i] = solver.eigenvalues()[i];
        }

        output = solver.eigenvectors();
    }

    auto canorg(const S_matrix& S) -> X_matrix
    {
        uint32_t size = S.rows();
        Eigen::MatrixXd eigenvectors;
        std::vector<double> eigenvalues;
        compute_eigenvalues(S, eigenvectors, eigenvalues);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(size, size);
        for (uint32_t i = 0; i < size; i++)
        {
            D(i, i) = 1.0 / std::sqrt(eigenvalues[i]);
        }

        return eigenvectors * (D * eigenvectors.transpose());
    }

    constexpr auto fact_ratio2(int64_t a, int64_t b) -> double { return fact(a) / fact(b) / fact(a - 2 * b); }

    constexpr auto b0(int i, int r, double g) -> double { return fact_ratio2(i, r) * pow_i(4 * g, r - i); }

    constexpr auto force_b(int64_t i, int64_t l1, int64_t l2, double p, double a, double b, int64_t r, double g) -> double
    {
        return binomial_prefactor(i, l1, l2, p - a, p - b) * b0(i, r, g);
    }

    constexpr auto b_term(int64_t i1, int64_t i2, int64_t r1, int64_t r2, int64_t u, int64_t l1, int64_t l2, int64_t l3, int64_t l4, double px,
                          double ax, double bx, double qx, double cx, double dx, double gamma1, double gamma2, double delta) -> double
    {
        return force_b(i1, l1, l2, px, ax, bx, r1, gamma1) * pow_i(-1, i2) * force_b(i2, l3, l4, qx, cx, dx, r2, gamma2) * pow_i(-1, u) *
               fact_ratio2(i1 + i2 - 2 * (r1 + r2), u) * pow_i(qx - px, i1 + i2 - 2 * (r1 + r2) - 2 * u) / pow_i(delta, i1 + i2 - 2 * (r1 + r2) - u);
    }

    constexpr auto b_array(int64_t l1, int64_t l2, int64_t l3, int64_t l4, double p, double a, double b, double q, double c, double d, double g1,
                           double g2, double delta) -> std::vector<double>
    {
        int imax = l1 + l2 + l3 + l4 + 1;
        std::vector<double> arr(imax, 0);
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

    constexpr auto a_term(int64_t i, int64_t r, int64_t u, int64_t l1, int64_t l2, double pax, double pbx, double cpx, double gamma) -> double
    {
        return pow_i(-1, i) * binomial_prefactor(i, l1, l2, pax, pbx) * pow_i(-1, u) * fact(i) * pow_i(cpx, i - 2 * r - 2 * u) *
               pow_i(0.25 / gamma, r + u) / fact(r) / fact(u) / fact(i - 2 * r - 2 * u);
    }

    constexpr auto a_array(int64_t l1, int64_t l2, double pa, double pb, double cp, double g) -> std::vector<double>
    {
        int imax = l1 + l2 + 1;
        std::vector<double> arr(imax, 0);
        for (int i = 0; i < imax; i++)
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

    constexpr auto overlap_1d(int64_t l1, int64_t l2, double x1, double x2, double gamma) -> double
    {
        double sum = 0;
        for (int i = 0; i < (1 + floor(0.5 * (l1 + l2))); i++)
        {
            sum += binomial_prefactor(2 * i, l1, l2, x1, x2) * fact2(2 * i - 1) / pow_i(2 * gamma, i);
        }
        return sum;
    }

    auto overlap(double alpha1, int l1, int m1, int n1, const glm::dvec3& a, double alpha2, int l2, int m2, int n2, const glm::dvec3& b) -> double
    {
        double rab2 = glm::distance2(a, b);
        double gamma = alpha1 + alpha2;
        glm::dvec3 p = gaussian_product_center(alpha1, a, alpha2, b);
        double pre = pow(M_PI / gamma, 1.5) * exp(-alpha1 * alpha2 * rab2 / gamma);
        double wx = overlap_1d(l1, l2, p.x - a.x, p.x - b.x, gamma);
        double wy = overlap_1d(m1, m2, p.y - a.y, p.y - b.y, gamma);
        double wz = overlap_1d(n1, n2, p.z - a.z, p.z - b.z, gamma);
        return pre * wx * wy * wz;
    }

    auto gto_overlap(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
    {
        return overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), pos1, gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(),
                       pos2);
    }

    auto cgf_overlap(const contracted_gaussian_functions& cgf1, const contracted_gaussian_functions& cgf2) -> double
    {
        uint32_t i = cgf1.orbs().size();
        uint32_t j = cgf2.orbs().size();
        double sum = 0;
        for (uint32_t k = 0; k < i; k++)
        {
            for (uint32_t l = 0; l < j; l++)
            {
                sum +=
                    cgf1.orbs()[k].factor() * cgf2.orbs()[l].factor() * gto_overlap(cgf1.orbs()[k], cgf2.orbs()[l], cgf1.get_pos(), cgf2.get_pos());
            }
        }
        return sum;
    }

    auto gto_kinetic(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
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
        uint32_t i = cgf1.orbs().size();
        uint32_t j = cgf2.orbs().size();
        double sum = 0;
        for (uint32_t k = 0; k < i; k++)
        {
            for (uint32_t l = 0; l < j; l++)
            {
                sum +=
                    cgf1.orbs()[k].factor() * cgf2.orbs()[l].factor() * gto_kinetic(cgf1.orbs()[k], cgf2.orbs()[l], cgf1.get_pos(), cgf2.get_pos());
            }
        }
        return sum;
    }

    auto nuclear(const glm::dvec3& a, double norm1, int l1, int m1, int n1, double alpha1, const glm::dvec3& b, double norm2, int l2, int m2, int n2,
                 double alpha2, const glm::dvec3& c) -> double
    {
        double gamma = alpha1 + alpha2;
        glm::dvec3 p = gaussian_product_center(alpha1, a, alpha2, b);
        double rab2 = glm::distance2(a, b);
        double rcp2 = glm::distance2(c, p);
        std::vector<double> ax = a_array(l1, l2, p.x - a.x, p.x - b.x, p.x - c.x, gamma);
        std::vector<double> ay = a_array(m1, m2, p.y - a.y, p.y - b.y, p.y - c.y, gamma);
        std::vector<double> az = a_array(n1, n2, p.z - a.z, p.z - b.z, p.z - c.z, gamma);
        double sum = 0.0;
        for (int i = 0; i <= l1 + l2; i++)
        {
            for (int j = 0; j <= m1 + m2; j++)
            {
                for (int k = 0; k <= n1 + n2; k++)
                {
                    sum += ax[i] * ay[j] * az[k] * Fgamma(i + j + k, rcp2 * gamma);
                }
            }
        }
        return -norm1 * norm2 * 2 * M_PI / gamma * exp(-alpha1 * alpha2 * rab2 / gamma) * sum;
    }

    auto gto_nuclear(const gto_data& gto1, const gto_data& gto2, const glm::dvec3& c, const glm::dvec3& pos1, const glm::dvec3& pos2) -> double
    {
        return nuclear(pos1, gto1.get_norm(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(), pos2, gto2.get_norm(), gto2.get_l(),
                       gto2.get_m(), gto2.get_n(), gto2.get_alpha(), c);
    }

    auto cgf_nuclear(contracted_gaussian_functions& cgf1, contracted_gaussian_functions& cgf2, const atom& a) -> double
    {
        uint32_t i = cgf1.orbs().size();
        uint32_t j = cgf2.orbs().size();
        double sum = 0;
#pragma omp parallel for collapse(2) reduction(+ : sum)
        for (uint32_t k = 0; k < i; k++)
        {
            for (uint32_t l = 0; l < j; l++)
            {
                sum += cgf1.orbs()[k].get_coeff() * cgf2.orbs()[l].get_coeff() *
                       gto_nuclear(cgf1.orbs()[k], cgf2.orbs()[l], a.pos(), cgf1.get_pos(), cgf2.get_pos());
            }
        }
        return sum * a.proton_count();
    }

    auto coulomb_repulsion(const glm::dvec3& a, double norma, int64_t la, int64_t ma, int64_t na, double alphaa, const glm::dvec3& b, double normb,
                           int64_t lb, int64_t mb, int64_t nb, double alphab, const glm::dvec3& c, double normc, int64_t lc, int64_t mc, int64_t nc,
                           double alphac, const glm::dvec3& d, double normd, int64_t ld, int64_t md, int64_t nd, double alphad) -> double
    {
        double rab2 = glm::distance2(a, b);
        double rcd2 = glm::distance2(c, d);
        glm::dvec3 p = gaussian_product_center(alphaa, a, alphab, b);
        glm::dvec3 q = gaussian_product_center(alphac, c, alphad, d);
        double rpq2 = glm::distance2(p, q);
        double gamma1 = alphaa + alphab;
        double gamma2 = alphac + alphad;
        double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);
        std::vector<double> bx = b_array(la, lb, lc, ld, p.x, a.x, b.x, q.x, c.x, d.x, gamma1, gamma2, delta);
        std::vector<double> by = b_array(ma, mb, mc, md, p.y, a.y, b.y, q.y, c.y, d.y, gamma1, gamma2, delta);
        std::vector<double> bz = b_array(na, nb, nc, nd, p.z, a.z, b.z, q.z, c.z, d.z, gamma1, gamma2, delta);
        double sum = 0.0;
        for (int i = 0; i <= (la + lb + lc + ld); i++)
        {
            for (int j = 0; j <= (ma + mb + mc + md); j++)
            {
                for (int k = 0; k <= (na + nb + nc + nd); k++)
                {
                    sum += bx[i] * by[j] * bz[k] * Fgamma(i + j + k, 0.25 * rpq2 / delta);
                }
            }
        }

        static constexpr auto SQRT_PI_FIFTH = 17.493418327624862;

        return 2 * SQRT_PI_FIFTH / (gamma1 * gamma2 * sqrt(gamma1 + gamma2)) * exp(-alphaa * alphab * rab2 / gamma1) *
               exp(-alphac * alphad * rcd2 / gamma2) * sum * norma * normb * normc * normd;
    }

    auto gto_repulsion(const gto_data& gto1, const gto_data& gto2, const gto_data& gto3, const gto_data& gto4, const glm::dvec3& pos1,
                       const glm::dvec3& pos2, const glm::dvec3& pos3, const glm::dvec3& pos4) -> double
    {
        return coulomb_repulsion(pos1, gto1.get_norm(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(), pos2, gto2.get_norm(),
                                 gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), pos3, gto3.get_norm(), gto3.get_l(), gto3.get_m(),
                                 gto3.get_n(), gto3.get_alpha(), pos4, gto4.get_norm(), gto4.get_l(), gto4.get_m(), gto4.get_n(), gto4.get_alpha());
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

    auto te_index(uint32_t i, uint32_t j, uint32_t k, uint32_t l) -> uint32_t
    {
        if (i < j)
        {
            std::swap(i, j);
        }
        if (k < l)
        {
            std::swap(k, l);
        }
        uint32_t ij = i * (i + 1) / 2 + j;
        uint32_t kl = k * (k + 1) / 2 + l;
        if (ij < kl)
        {
            std::swap(ij, kl);
        }
        return ij * (ij + 1) / 2 + kl;
    }

} // namespace

class hartree_fock
{
public:
    std::shared_ptr<const molecule> mol;
    uint32_t orbital_count{};
    uint32_t atom_count{};
    uint32_t electron_count{};
    uint32_t step_count{};
    std::vector<std::string> orbital_desc_list;
    std::vector<contracted_gaussian_functions> orbitals;
    H_matrix H;
    S_matrix S;
    T_matrix T;
    P_matrix P;
    P_matrix P_new;
    G_matrix G;
    F_matrix F;
    F_matrix Fp;
    X_matrix X;
    X_matrix Xp;
    C_matrix C;
    Cc_matrix Cc;
    std::vector<V_matrix> V;
    std::vector<double> TE;
    std::vector<double> energies;
    std::vector<double> molorb_energy;
    double energy{};
    double nuclear_repulsion{};
    double alpha{0.5};

    hartree_fock() = default;

    void set_molecule(const std::shared_ptr<const molecule>& target)
    {
        mol = target;
        electron_count = 0;
        energy = 0;
        step_count = 0;
        uint32_t count = 0;
        for (uint32_t i = 0; i < mol->atoms_count(); i++)
        {
            electron_count += mol->get_atoms()[i].electron_count();
            for (uint32_t j = 0; j < mol->get_atoms()[i].orbital_count(); j++)
            {
                count++;
                orbital_desc_list.push_back(fmt::format("[{}] {}{}\t ({})", count, mol->get_atoms()[i].name(), i + 1, mol->get_atoms()[i][j].type()));
                orbitals.push_back(mol->get_atoms()[i][j]);
            }
        }
        electron_count -= mol->get_charge();
        orbital_count = count;
        atom_count = mol->atoms_count();
        S = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        H = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        T = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        V.resize(atom_count);
        P = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        Cc = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        C = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        X = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        Xp = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        G = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        for (auto& mat : V)
        {
            mat.resize(orbital_count, orbital_count);
        }
        TE.resize(te_index(orbital_count, orbital_count, orbital_count, orbital_count) + 1, -1.0);
    }

    auto run() -> hartree_fock_result
    {
        setup();
        uint32_t iter = iterate();
        auto atoms = mol->get_atoms() | std::views::transform([](const atom& ato) { return std::make_pair(ato.proton_count(), ato.pos()); });
        return {.iterations = iter,
                .mo_energies = std::move(molorb_energy),
                .coefficients = std::move(C),
                .orbitals = std::move(orbitals),
                .atoms = std::vector<std::pair<uint32_t, glm::dvec3>>(atoms.begin(), atoms.end()),
                .electron_count = electron_count,
                .homo_index = (electron_count + 1) / 2 - 1,
                .basis_name = mol->get_basis_name()};
    }

    void step()
    {
        step_count++;
        uint32_t index1 = 0;
        uint32_t index2 = 0;
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                G(i, j) = 0.;
                for (uint32_t k = 0; k < orbital_count; k++)
                {
                    for (uint32_t l = 0; l < orbital_count; l++)
                    {
                        index1 = te_index(i, j, l, k);
                        index2 = te_index(i, k, l, j);
                        G(i, j) += P(k, l) * (TE[index1] - 0.5 * TE[index2]);
                    }
                }
            }
        }
        F = H + G;
        Fp = Xp * F * X;
        compute_eigenvalues(Fp, Cc, molorb_energy);
        energy = calc_energy();
        C = X * Cc;
        P_new = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                for (uint32_t k = 0; k < electron_count / 2; k++)
                {
                    P_new(i, j) += 2.0 * C(i, k) * C(j, k);
                }
            }
        }
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                P(i, j) = (1.0 - alpha) * P_new(i, j) + alpha * P(i, j);
            }
        }
    }

    void setup()
    {

#pragma omp parallel for
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                S(i, j) = cgf_overlap(orbitals[i], orbitals[j]);
            }
        }

#pragma omp parallel for
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                T(i, j) = cgf_kinetic(orbitals[i], orbitals[j]);
                H(i, j) = T(i, j);
            }
        }

#pragma omp parallel for
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                for (uint32_t k = 0; k < atom_count; k++)
                {
                    V[k](i, j) = cgf_nuclear(orbitals[i], orbitals[j], mol->get_atoms()[k]);
                    H(i, j) += V[k](i, j);
                }
            }
        }

        std::vector<std::array<uint32_t, 5>> te_jobs;
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j <= i; j++)
            {
                uint32_t ij = i * (i + 1) / 2 + j;
                for (uint32_t k = 0; k < orbital_count; k++)
                {
                    for (uint32_t l = 0; l <= k; l++)
                    {
                        uint32_t kl = k * (k + 1) / 2 + l;
                        if (ij <= kl)
                        {
                            uint32_t index = te_index(i, j, k, l);
                            te_jobs.push_back({index, i, j, k, l});
                        }
                    }
                }
            }
        }
#pragma omp parallel for schedule(dynamic)
        for (const auto& te_job : te_jobs)
        {
            TE[te_job[0]] = cgf_repulsion(orbitals[te_job[1]], orbitals[te_job[2]], orbitals[te_job[3]], orbitals[te_job[4]]);
        }

        X = canorg(S);
        Xp = X.transpose();
        nuclear_repulsion = calculate_nuclear_repulsion();
    }

    auto iterate() -> uint32_t
    {
        double energy_diff = 1;
        double oldenergy = 1000;
        uint32_t iter = 0;

        while (energy_diff > 1e-9)
        {
            iter++;
            step();
            energy_diff = std::abs(energy - oldenergy);
            oldenergy = energy;
            energies.push_back(energy);

            if (iter > 1000)
            {
                std::cout << "Too may iteration steps.. terminating and outputting results." << std::endl;
                break;
            }
        }

        return iter;
    }

    auto calc_energy() -> double
    {
        double energy = 0;
        Eigen::MatrixXd mat = (H * F);
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                energy += P(j, i) * mat(i, j);
            }
        }
        return energy * 0.5 + nuclear_repulsion;
    }

    auto calculate_nuclear_repulsion() -> double
    {
        double repulsion = 0;
        for (uint32_t i = 0; i < atom_count; i++)
        {
            for (uint32_t j = i + 1; j < atom_count; j++)
            {
                repulsion += mol->get_atoms()[i].proton_count() * mol->get_atoms()[j].proton_count() /
                             glm::length(mol->get_atoms()[i].pos() - mol->get_atoms()[j].pos());
            }
        }
        return repulsion;
    }
};

auto solve(const std::shared_ptr<molecule>& molecule) -> hartree_fock_result
{
    hartree_fock solver;
    solver.set_molecule(molecule);
    return solver.run();
}

void write_result(const hartree_fock_result& result, nlohmann::json& out_json)
{
    size_t mat_rows = result.coefficients.rows();
    std::vector<double> coefficient_list(mat_rows * mat_rows);

    for (size_t i = 0; i < mat_rows; i++)
    {
        for (size_t j = 0; j < mat_rows; j++)
        {
            coefficient_list[i * mat_rows + j] = result.coefficients(i, j);
        }
    }

    nlohmann::json atoms(nlohmann::json::array_t{});

    for (const auto& inst : result.atoms)
    {
        atoms.push_back({
            inst.first,
            inst.second.x,
            inst.second.y,
            inst.second.z,
        });
    }

    out_json = {{"iterations", result.iterations},
                {"mo_energies", result.mo_energies},
                {"coefficient_n", mat_rows},
                {"coefficients", coefficient_list},
                {"atoms", atoms},
                {"electron_count", result.electron_count},
                {"homo_index", result.homo_index},
                {"basis_name", result.basis_name},
                {"type", "mo_output"}};
}

void read_result(hartree_fock_result& result, const nlohmann::json& in_json)
{
    result = {};

    result.iterations = in_json.at("iterations");
    result.mo_energies = in_json.at("mo_energies").get<std::vector<double>>();
    result.electron_count = in_json.at("electron_count");
    result.homo_index = in_json.at("homo_index");
    result.basis_name = in_json.at("basis_name");

    size_t mat_rows = in_json.at("coefficient_n");

    result.coefficients.resize(mat_rows, mat_rows);
    for (size_t i = 0; i < mat_rows; i++)
    {
        for (size_t j = 0; j < mat_rows; j++)
        {
            result.coefficients(i, j) = in_json.at("coefficients")[i * mat_rows + j];
        }
    }

    auto basis = basis_manager::get_instance().get_basis(in_json.at("basis_name"));
    for (const auto& [_, value] : in_json.at("atoms").items())
    {
        auto atomic_number = value[0].get<uint32_t>();
        glm::dvec3 pos = {value[1].get<double>(), value[2].get<double>(), value[3].get<double>()};
        result.atoms.emplace_back(atomic_number, pos);

        for (const auto& atomic_no : basis->atom_data(atomic_number))
        {
            result.orbitals.emplace_back(pos, atomic_no);
        }
    }
}
