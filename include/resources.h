#pragma once

#include "singleton.h"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

class resource_manager : public singleton<resource_manager>
{
    std::unordered_map<std::string, std::vector<uint8_t>> resources;
protected:
    resource_manager();
public:

    inline auto get_resource(const std::string& name) const -> const auto& { return resources.at(name); }
};
