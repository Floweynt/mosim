#include "hf/chem.h"

auto molecule::read_from_file(const std::string& file) -> std::shared_ptr<molecule>
{
    std::ifstream ifs(file);
    if (!ifs)
    {
        throw std::runtime_error("failed to open file: " + file);
    }
    nlohmann::json config = nlohmann::json::parse(ifs);
    auto instance = std::make_shared<molecule>(basis_manager::get_instance().get_basis(config.value("basis", "sto-6g")));

    double units = 1;
    std::string unit_enumerator = config.value("units", "bohr");
    if (unit_enumerator == "bohr")
    {
        units = 1;
    }
    else if (unit_enumerator == "angstrom")
    {
        units = 1.88973;
    }
    else if (unit_enumerator == "nm")
    {
        units = 18.8973;
    }
    else
    {
        throw std::runtime_error("bad unit type: " + unit_enumerator);
    }

    instance->charge = config.value("charge", uint32_t(0));

    for (const auto& atom_spec : config.at("atoms").items())
    {
        std::string element = atom_spec.value().at(0).get<std::string>();
        instance->add_atom(element, atom_spec.value().at(1).get<double>() * units, atom_spec.value().at(2).get<double>() * units,
                           atom_spec.value().at(3).get<double>() * units);
    }
    return instance;
}

