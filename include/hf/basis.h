#pragma once
#include "math_util.h"
#include "singleton.h"
#include "util.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <glm/fwd.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class gto_data
{
    double coeff;
    uint32_t l, m, n;
    double alpha;
    double norm;
    double coeff_norm_factor;

    [[nodiscard]] auto calc_norm() const -> double
    {
        double nom = std::pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) * std::pow(alpha, (l + m + n) + 3.0 / 2.0);
        double denom =
            (double)fact2(2 * int64_t(l) - 1) * (double)fact2(2 * int64_t(m) - 1) * (double)fact2(2 * int64_t(n) - 1) * std::pow(M_PI, 3.0 / 2.0);
        return std::sqrt(nom / denom);
    }

public:
    gto_data(double coeff, uint32_t l, uint32_t m, uint32_t n, double alpha)
        : coeff(coeff), l(l), m(m), n(n), alpha(alpha), norm(calc_norm()), coeff_norm_factor(coeff * norm)
    {
    }

    [[nodiscard]] constexpr auto get_coeff() const { return coeff; }
    [[nodiscard]] constexpr auto get_l() const { return l; }
    [[nodiscard]] constexpr auto get_m() const { return m; }
    [[nodiscard]] constexpr auto get_n() const { return n; }
    [[nodiscard]] constexpr auto get_alpha() const { return alpha; }
    [[nodiscard]] constexpr auto get_norm() const { return norm; }
    [[nodiscard]] constexpr auto factor() const { return coeff_norm_factor; }
};

class cgf_data
{
    std::string type;
    std::vector<gto_data> data;

public:
    cgf_data(std::string type, std::vector<gto_data> data) : type(std::move(type)), data(std::move(data)) {}

    [[nodiscard]] constexpr auto get_type() const -> const auto& { return type; }
    [[nodiscard]] constexpr auto orbs() const -> const auto& { return data; }
};

class basis_set
{
    // atoms<orbitals<N>>
    std::unordered_map<uint32_t, std::vector<cgf_data>> dataset;

    bool initalized{};

    basis_set() = default;
    basis_set(const std::string& file) { load(file); }

    void init_json(const nlohmann::json& config);

    void load(const std::string& file);
    void load_from_builtin(const std::string& file);

public:
    [[nodiscard]] constexpr auto is_initalized() const { return initalized; }
    auto atom_data(uint32_t atomic_number) const -> const auto&
    {
        if (!initalized)
        {
            throw std::runtime_error("dataset not loaded");
        }
        return dataset.at(atomic_number);
    }

    static auto load_basis(const std::string& path) { return std::make_shared<basis_set>(basis_set(path)); }
    static auto load_basis_builtin(const std::string& path)
    {
        auto value = std::make_shared<basis_set>(basis_set());
        value->load_from_builtin(path);
        return value;
    }
};

class basis_manager : public singleton<basis_manager>
{
    std::unordered_map<std::string, std::shared_ptr<basis_set>> basis_list;

public:
    void register_basis(const std::string& name) { basis_list[name] = basis_set::load_basis_builtin("assets/basis/" + name + ".json"); }
    inline auto get_basis(const std::string& name) const -> const auto& { return basis_list.at(name); }
};
