#include "hf/chem.h"
#include <fmt/core.h>

inline static constexpr auto ANGSTROM_TO_BOHR = 1.88973;
inline static constexpr auto NM_TO_BOHR = 18.8973;

namespace
{
    inline constexpr std::pair<const char*, double> UNIT_LUT[] = {
        {"bohr", 1}, {"angstrom", ANGSTROM_TO_BOHR}, {"A", ANGSTROM_TO_BOHR}, {"nanometer", NM_TO_BOHR}, {"nm", NM_TO_BOHR}};
    auto parse_unit(const std::string_view& name) -> double
    {
        for (const auto& unit_desc : UNIT_LUT)
        {
            if (name == unit_desc.first)
            {
                return unit_desc.second;
            }
        }

        throw std::runtime_error(fmt::format("unknown unit: {}", name));
    }
} // namespace

auto molecule::read_from_file(const nlohmann::json& config) -> std::shared_ptr<molecule>
{
    std::string basis_name;
    auto instance = std::make_shared<molecule>(basis_manager::get_instance().get_basis(basis_name = config.value("basis", "sto-6g")));
    instance->basis_name = basis_name;
    double units = parse_unit(config.value("units", "bohr"));

    instance->charge = config.value("charge", uint32_t(0));

    for (const auto& atom_data : config.at("atoms").items())
    {
        instance->add_atom(atom_data.value().at(0).get<std::string>(), atom_data.value().at(1).get<double>() * units,
                           atom_data.value().at(2).get<double>() * units, atom_data.value().at(3).get<double>() * units);
    }
    return instance;
}

