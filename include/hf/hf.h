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
};

auto solve(const std::shared_ptr<molecule>& molecule) -> hartree_fock_result;
