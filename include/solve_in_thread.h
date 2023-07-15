#pragma once

#include "hf/hf.h"
#include <chrono>
#include <fstream>
#include <future>
#include <nlohmann/json.hpp>

inline auto run_hf(const std::string& path) -> std::future<hartree_fock_result>
{
    return std::async(std::launch::async, [path]() {
        std::ifstream ifs(path);
        auto json = nlohmann::json::parse(ifs);
        hartree_fock_result result;

        if (json.value("type", "") == "mo_output")
        {
            read_result(result, json);
        }
        else
        {
            auto mol = molecule::read_from_file(json);
            result = solve(mol);
        }

        return result;
    });
}

template <typename T>
auto is_ready(const std::future<T>& future) -> bool
{
    return future.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}
