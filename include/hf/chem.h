// cSpell:ignore nlohmann
#pragma once

#include "basis.h"
#include "math_util.h"
#include <array>
#include <cstdint>
#include <fmt/core.h>
#include <fstream>
#include <glm/glm.hpp>
#include <iomanip>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <vector>

class contracted_gaussian_functions
{
    const cgf_data* data;
    glm::dvec3 pos;

public:
    constexpr contracted_gaussian_functions(const glm::dvec3& pos, const cgf_data& data) : data(&data), pos(pos) {}

    [[nodiscard]] constexpr auto get_data() const -> const auto& { return *data; }
    [[nodiscard]] constexpr auto get_pos() const -> const auto& { return pos; }

    [[nodiscard]] constexpr auto type() const -> const auto& { return get_data().get_type(); }
    [[nodiscard]] constexpr auto orbs() const -> const auto& { return get_data().orbs(); }
};

inline static constexpr std::array ELEMENT_NAMES = {"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
                                                    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
                                                    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"};

inline static constexpr std::array ELEMENT_COLORS = {0xffffff, 0xd9ffff, 0xcc80ff, 0xc2ff00, 0xffb5b5, 0x909090, 0x3050f8, 0xff0d0d, 0x90e050,
                                                     0xb3e3f5, 0xab5cf2, 0x8aff00, 0xbfa6a6, 0xf0c8a0, 0xff8000, 0xffff30, 0x1ff01f, 0x80d1e3,
                                                     0x8f40d4, 0x3dff00, 0xe6e6e6, 0xbfc2c7, 0xa6a6ab, 0x8a99c7, 0x9c7ac7, 0xe06633, 0xf090a0,
                                                     0x50d050, 0xc88033, 0x7d80b0, 0xc28f8f, 0x668f8f, 0xbd80e3, 0xffa100, 0xa62929, 0x5cb8d1};

class atom
{
private:
    uint32_t atomic_number;
    glm::dvec3 r{};
    std::string element;
    std::vector<contracted_gaussian_functions> wavefunctions;
    uint32_t electrons;

public:
    atom(const std::string& element_name, double x, double y, double z, const basis_set& basis)
        : atomic_number(e2z(element_name)), r(x, y, z), element(element_name), electrons(atomic_number)
    {
        add_wavefunctions(basis);
    }

     atom(uint32_t atomic_number, double x, double y, double z, const basis_set& basis)
        : atomic_number(atomic_number), r(x, y, z), electrons(atomic_number)
    {
        add_wavefunctions(basis);
    }

private:
    [[nodiscard]] constexpr auto x() const -> double { return r.x; };
    [[nodiscard]] constexpr auto y() const -> double { return r.y; };
    [[nodiscard]] constexpr auto z() const -> double { return r.z; };

    static constexpr auto e2z(const std::string& symbol) -> uint32_t
    {
        size_t index = 0;
        for (const auto* name : ELEMENT_NAMES)
        {
            if (name == symbol)
            {
                return index + 1;
            }

            index++;
        }

        throw std::runtime_error("Unknown element: " + symbol);
    }

    static constexpr auto z2e(uint32_t atomic_number) -> std::string_view
    {
        atomic_number -= 1;
        if (atomic_number >= ELEMENT_NAMES.size())
        {
            throw std::runtime_error(fmt::format("Unknown element: {}", atomic_number));
        }
        return ELEMENT_NAMES[atomic_number];
    }

    void add_wavefunctions(const basis_set& basis)
    {
        const auto& atom_data = basis.atom_data(atomic_number);

        for (const auto& orbital : atom_data)
        {
            wavefunctions.emplace_back(r, orbital);
        }
    }

public:
    [[nodiscard]] constexpr auto name() const -> std::string_view { return z2e(atomic_number); }
    [[nodiscard]] constexpr auto pos() const -> const auto& { return r; }
    [[nodiscard]] constexpr auto orbital_count() const -> uint32_t { return wavefunctions.size(); }
    constexpr auto operator[](uint32_t index) const -> const contracted_gaussian_functions& { return wavefunctions[index]; }
    [[nodiscard]] constexpr auto electron_count() const -> uint32_t { return electrons; }
    [[nodiscard]] constexpr auto proton_count() const -> uint32_t { return atomic_number; }
};

class molecule
{
private:
    std::vector<atom> atoms;
    std::shared_ptr<const basis_set> basis;
    uint32_t charge{};
    std::string basis_name;

public:
    molecule(const std::shared_ptr<const basis_set>& basis) : basis(basis){};

    void add_atom(const atom& at) { atoms.push_back(at); }
    void add_atom(const std::string& symbol, double x, double y, double z) { atoms.emplace_back(symbol, x, y, z, *basis); }
    void add_atom(const uint32_t atomic_number, double x, double y, double z) { atoms.emplace_back(atomic_number, x, y, z, *basis); }

    [[nodiscard]] constexpr auto atoms_count() const { return atoms.size(); };
    [[nodiscard]] constexpr auto get_atoms() const -> const auto& { return atoms; }
    [[nodiscard]] constexpr auto get_atoms() -> auto& { return atoms; }
    [[nodiscard]] constexpr auto get_charge() const -> const auto& { return charge; }
    [[nodiscard]] constexpr auto get_basis_name() const -> const auto& { return basis_name; }

    static auto read_from_file(const nlohmann::json& file) -> std::shared_ptr<molecule>;
};

