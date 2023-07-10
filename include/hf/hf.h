#pragma once

#include "chem.h"
#include <eigen3/Eigen/Eigen>

struct hartree_fock_result
{
    size_t iterations;
    std::vector<double> mo_energies;
    Eigen::MatrixXd coefficients;
    std::vector<contracted_gaussian_functions> orbitals;
    std::vector<std::pair<uint32_t, glm::dvec3>> atoms;
    size_t electron_count;
    size_t homo_index;
    std::string basis_name;
};

auto solve(const std::shared_ptr<molecule>& molecule) -> hartree_fock_result;
void write_result(const hartree_fock_result& result, nlohmann::json& out_json);
void read_result(hartree_fock_result& result, const nlohmann::json& in_json);
