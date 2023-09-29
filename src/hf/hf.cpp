#include "hf/hf.h"
#include "Eigen/src/Core/Matrix.h"
#include "hf/basis.h"
#include "hf/chem.h"
#include "hf/hf_setup.h"
#include <fmt/core.h>
#include <iostream>
#include <mutex>
#include <ranges>
#include <unordered_set>

using matrix = Eigen::MatrixXd;
using vect = Eigen::VectorXd;

namespace
{
    auto diagonalize(const matrix& mat)
    {
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
        return std::make_pair(solver.eigenvalues(), solver.eigenvectors());
    }

    auto canorg(const matrix& S) -> matrix
    {
        auto [eigenvalues, eigenvectors] = diagonalize(S);
        auto D = eigenvalues.cwiseSqrt().cwiseInverse().asDiagonal();
        return eigenvectors * (D * eigenvectors.transpose());
    }
} // namespace

namespace detail
{
    struct pairhash
    {
    public:
        template <typename T, typename U>
        std::size_t operator()(const std::pair<T, U>& x) const
        {
            return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
        }
    };
    static std::unordered_set<std::pair<int64_t, double>, pairhash> tab;
    static std::atomic<size_t> count = 0;

    struct spinlock
    {
        std::atomic<bool> lock_ = {0};

        void lock() noexcept
        {
            for (;;)
            {
                // Optimistically assume the lock is free on the first try
                if (!lock_.exchange(true, std::memory_order_acquire))
                {
                    return;
                }
                // Wait for lock to be released without generating cache misses
                while (lock_.load(std::memory_order_relaxed))
                {
                    // Issue X86 PAUSE or ARM YIELD instruction to reduce contention between
                    // hyper-threads
                    __builtin_ia32_pause();
                }
            }
        }

        bool try_lock() noexcept
        {
            // First do a relaxed load to check if lock is free in order to prevent
            // unnecessary cache misses if someone does while(!try_lock())
            return !lock_.load(std::memory_order_relaxed) && !lock_.exchange(true, std::memory_order_acquire);
        }

        void unlock() noexcept { lock_.store(false, std::memory_order_release); }
    };

    spinlock lock;
    void report(int64_t m, double x)
    {
        count++;
        std::lock_guard g(lock);
        tab.emplace(m, x);
    }
} // namespace detail

class hartree_fock
{
public:
    std::shared_ptr<const molecule> mol;
    uint32_t orbital_count{};
    uint32_t atom_count{};
    uint32_t electron_count{};
    std::vector<contracted_gaussian_functions> orbitals;

    matrix hamiltonian;
    matrix overlap;
    matrix kinetic;
    std::vector<double> two_electrons;
    std::vector<matrix> nuclear;

    matrix P;
    matrix P_new;
    matrix G;
    matrix F;
    matrix Fp;
    matrix X;
    matrix Xp;
    matrix C;
    matrix Cc;

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
        uint32_t count = 0;
        for (uint32_t i = 0; i < mol->atoms_count(); i++)
        {
            electron_count += mol->get_atoms()[i].electron_count();
            for (uint32_t j = 0; j < mol->get_atoms()[i].orbital_count(); j++)
            {
                count++;
                orbitals.push_back(mol->get_atoms()[i][j]);
            }
        }
        electron_count -= mol->get_charge();
        orbital_count = count;
        atom_count = mol->atoms_count();
        overlap = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        hamiltonian = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        kinetic = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        nuclear.resize(atom_count);
        P = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        Cc = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        C = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        X = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        Xp = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        G = Eigen::MatrixXd::Zero(orbital_count, orbital_count);
        for (auto& mat : nuclear)
        {
            mat.resize(orbital_count, orbital_count);
        }
        two_electrons.resize(hf::te_index(orbital_count, orbital_count, orbital_count, orbital_count) + 1, -1.0);
    }

    auto run() -> hartree_fock_result
    {
        auto start = std::chrono::high_resolution_clock::now();
        setup();
        uint32_t iter = iterate();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(stop - start);

        auto atoms = mol->get_atoms() | std::views::transform([](const atom& ato) { return std::make_pair(ato.proton_count(), ato.pos()); });
        std::vector<std::string> ao_names;

        for (uint32_t i = 0; i < mol->atoms_count(); i++)
        {
            for (uint32_t j = 0; j < mol->get_atoms()[i].orbital_count(); j++)
            {
                ao_names.emplace_back(fmt::format("{}{}-{}", mol->get_atoms()[i].name(), i, mol->get_atoms()[i][j].type()));
            }
        }

        std::cout << detail::count << " " <<  detail::tab.size() << '\n'; 
        return {.iterations = iter,
            .mo_energies = std::move(molorb_energy),
            .coefficients = std::move(C),
            .orbitals = std::move(orbitals),
            .atoms = std::vector<std::pair<uint32_t, glm::dvec3>>(atoms.begin(), atoms.end()),
            .electron_count = electron_count,
            .homo_index = (electron_count + 1) / 2 - 1,
            .basis_name = mol->get_basis_name(),
            .solve_time_microseconds = (size_t)duration.count(),
            .orbital_names = ao_names};
    }

    void step()
    {
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                G(i, j) = 0.;
                for (uint32_t k = 0; k < orbital_count; k++)
                {
                    for (uint32_t l = 0; l < orbital_count; l++)
                    {
                        G(i, j) += P(k, l) * (two_electrons[hf::te_index(i, j, l, k)] - 0.5 * two_electrons[hf::te_index(i, k, l, j)]);
                    }
                }
            }
        }
        F = hamiltonian + G;

        Fp = Xp * F * X;
        auto [eigen_val, eigen_vec] = diagonalize(Fp);
        Cc = eigen_vec;
        molorb_energy = std::vector(eigen_val.begin(), eigen_val.end());

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
        // calculate overlap matrix
#pragma omp parallel for
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                overlap(i, j) = hf::cgf_overlap(orbitals[i], orbitals[j]);
            }
        }

        // calculate kinetic matrix
#pragma omp parallel for
        for (uint32_t i = 0; i < orbital_count; i++)
        {
            for (uint32_t j = 0; j < orbital_count; j++)
            {
                kinetic(i, j) = hf::cgf_kinetic(orbitals[i], orbitals[j]);
            }
        }

        hamiltonian = kinetic;

        // calculate nuclear matrix
#pragma omp parallel for
        for (uint32_t k = 0; k < atom_count; k++)
        {
            for (uint32_t i = 0; i < orbital_count; i++)
            {
                for (uint32_t j = 0; j < orbital_count; j++)
                {
                    nuclear[k](i, j) = hf::cgf_nuclear(orbitals[i], orbitals[j], mol->get_atoms()[k]);
                }
            }
            hamiltonian += nuclear[k];
        }

        // calculate two-electron integrals
        std::vector<std::array<uint32_t, 5>> electron_jobs;
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
                            uint32_t index = hf::te_index(i, j, k, l);
                            electron_jobs.push_back({index, i, j, k, l});
                        }
                    }
                }
            }
        }
#pragma omp parallel for schedule(dynamic)
        for (const auto& te_job : electron_jobs)
        {
            two_electrons[te_job[0]] = hf::cgf_repulsion(orbitals[te_job[1]], orbitals[te_job[2]], orbitals[te_job[3]], orbitals[te_job[4]]);
        }

        X = canorg(overlap);
        Xp = X.transpose();
        nuclear_repulsion = calculate_nuclear_repulsion();
    }

    auto iterate() -> uint32_t
    {
        step();

        double energy_diff = 0;
        double oldenergy = energy;
        uint32_t iter = 0;

        while (true)
        {
            iter++;
            step();
            energy_diff = std::abs(energy - oldenergy);

            fmt::println("iter {}: delta={}, e={}", iter, energy_diff, energy);
            oldenergy = energy;
            energies.push_back(energy);

            if (energy_diff < 1e-9)
            {
                break;
            }

            if (iter > 10000)
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
        Eigen::MatrixXd mat = hamiltonian * F;
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

    out_json = {{"iterations", result.iterations}, {"mo_energies", result.mo_energies}, {"coefficient_n", mat_rows},
        {"coefficients", coefficient_list}, {"atoms", atoms}, {"electron_count", result.electron_count}, {"homo_index", result.homo_index},
        {"basis_name", result.basis_name}, {"orbital_names", result.orbital_names}, {"type", "mo_output"}};
}

void read_result(hartree_fock_result& result, const nlohmann::json& in_json)
{
    result = {};

    result.iterations = in_json.at("iterations");
    result.mo_energies = in_json.at("mo_energies").get<std::vector<double>>();
    result.electron_count = in_json.at("electron_count");
    result.homo_index = in_json.at("homo_index");
    result.basis_name = in_json.at("basis_name");
    result.orbital_names = in_json.at("orbital_names");

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
