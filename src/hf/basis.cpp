#include "hf/basis.h"
#include "resources.h"
#include <nlohmann/json_fwd.hpp>

inline static constexpr std::pair<std::string, glm::uvec3> ORBITALS_METADATA_S[] = {
    {"s", {0, 0, 0}},
};

inline static constexpr std::pair<std::string, glm::uvec3> ORBITALS_METADATA_P[] = {
    {"px", {1, 0, 0}},
    {"py", {0, 1, 0}},
    {"pz", {0, 0, 1}},
};

inline static constexpr std::pair<std::string, glm::uvec3> ORBITALS_METADATA_D[] = {
    {"dx2", {2, 0, 0}}, {"dxy", {1, 1, 0}}, {"dxz", {1, 0, 1}}, {"dy2", {0, 2, 0}}, {"dyz", {0, 1, 1}}, {"dz2", {0, 0, 2}},
};

inline static constexpr const std::pair<std::string, glm::uvec3>* ORBITALS_METADATA_LUT[] = {
    decay_val(ORBITALS_METADATA_S), decay_val(ORBITALS_METADATA_P), decay_val(ORBITALS_METADATA_D)};

inline static constexpr size_t ORBITALS_METADATA_LUT_SIZE[] = {
    std::extent<decltype(ORBITALS_METADATA_S)>::value,
    std::extent<decltype(ORBITALS_METADATA_P)>::value,
    std::extent<decltype(ORBITALS_METADATA_D)>::value,
};

void basis_set::init_json(const nlohmann::json& config)
{
    for (const auto& element_mo : config.at("elements").items())
    {
        std::vector<cgf_data> atom;
        auto shells = element_mo.value().at("electron_shells");
        size_t shell_no = 1;
        for (const auto& [_, shell] : shells.items())
        {
            std::vector<gto_data> shell_data;
            size_t coeff_count = shell.at("exponents").size();
            auto exponents = shell.at("exponents");

            for (size_t i = 0; i < shell.at("angular_momentum").size(); i++)
            {
                size_t angular = shell.at("angular_momentum").at(i);
                if (angular > 2)
                {
                    throw std::runtime_error("L > 2");
                }

                auto coefficients = shell.at("coefficients").at(i);

                const auto* metadata = ORBITALS_METADATA_LUT[angular];
                auto metadata_size = ORBITALS_METADATA_LUT_SIZE[angular];

                for (size_t metadata_index = 0; metadata_index < metadata_size; metadata_index++)
                {
                    const auto& orbital = metadata[metadata_index];
                    for (size_t j = 0; j < coeff_count; j++)
                    {
                        shell_data.emplace_back(std::stod(coefficients.at(j).get<std::string>()), orbital.second.x, orbital.second.y,
                                                orbital.second.z, std::stod(exponents.at(j).get<std::string>()));
                    }
                    atom.emplace_back(std::to_string(shell_no) + orbital.first, std::move(shell_data));
                }
            }

            shell_no++;
        }

        dataset.emplace(std::stoi(element_mo.key()), std::move(atom));
    }

    initalized = true;
}

void basis_set::load(const std::string& file)
{
    std::ifstream ifs(file);
    if (!ifs)
    {
        throw std::runtime_error("failed to open file: " + file);
    }

    nlohmann::json config = nlohmann::json::parse(ifs);
    init_json(config);
}

void basis_set::load_from_builtin(const std::string& file)
{
    const auto& resource = resource_manager::get_instance().get_resource(file);
    init_json(nlohmann::json::parse(std::string(resource.begin(), resource.end())));
}
