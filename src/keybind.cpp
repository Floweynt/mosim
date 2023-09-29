#include "keybind.h"
#include <nlohmann/json.hpp>
#include <string>

auto kb_from_json(const nlohmann::json& obj) -> keybind
{
    try
    {
        return keybind{.is_mouse = obj.at("is_mouse"), .key = obj.at("key"), .mods = obj.at("mods")};
    }
    catch (...)
    {
        return {};
    }
}

auto kb_to_json(const keybind& kb) -> nlohmann::json
{
    return {
        {"is_mouse", kb.is_mouse},
        {"key", kb.key},
        {"mods", kb.mods},
    };
}

void keybind_config::save(const std::string& path) const
{
    std::ofstream ofs(path, std::ios_base::trunc);
    if (!ofs)
    {
        throw std::runtime_error("failed to write to config");
        return;
    }

    ofs << nlohmann::json{
        {"goto_lumo", kb_to_json(goto_lumo)},
        {"goto_homo", kb_to_json(goto_homo)},
        {"goto_next_mo", kb_to_json(goto_next_mo)},
        {"goto_prev_mo", kb_to_json(goto_prev_mo)},
        {"goto_next_mo_alt", kb_to_json(goto_next_mo_alt)},
        {"goto_prev_mo_alt", kb_to_json(goto_prev_mo_alt)},
        {"open_config", kb_to_json(open_config)},
        {"toggle_render_mo", kb_to_json(toggle_render_mo)},
        {"toggle_render_atoms", kb_to_json(toggle_render_atoms)},
        {"toggle_render_box", kb_to_json(toggle_render_box)},
        {"increment_isolevel", kb_to_json(increment_isolevel)},
        {"decrement_isolevel", kb_to_json(decrement_isolevel)},
    };
}

void keybind_config::load(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs)
    {
        save(path);
        return;
    }

    try
    {
        nlohmann::json keybind = nlohmann::json::parse(ifs);
        goto_lumo = kb_from_json(keybind.at("goto_lumo"));
        goto_homo = kb_from_json(keybind.at("goto_homo"));
        goto_next_mo = kb_from_json(keybind.at("goto_next_mo"));
        goto_prev_mo = kb_from_json(keybind.at("goto_prev_mo"));
        goto_next_mo_alt = kb_from_json(keybind.at("goto_next_mo_alt"));
        goto_prev_mo_alt = kb_from_json(keybind.at("goto_prev_mo_alt"));
        open_config = kb_from_json(keybind.at("open_config"));
        toggle_render_mo = kb_from_json(keybind.at("toggle_render_mo"));
        toggle_render_atoms = kb_from_json(keybind.at("toggle_render_atoms"));
        toggle_render_box = kb_from_json(keybind.at("toggle_render_box"));
        increment_isolevel = kb_from_json(keybind.at("increment_isolevel"));
        decrement_isolevel = kb_from_json(keybind.at("decrement_isolevel"));
    }
    catch (...)
    {
        save(path);
    }
}

auto on_keybind(const keybind& keybind, int key, int mods, const std::function<void()>& handler, bool mouse) -> bool
{
    if (mouse != keybind.is_mouse)
    {
        return false;
    }

    if (key != keybind.key)
    {
        return false;
    }

    if ((mods & keybind.mods) != keybind.mods)
    {
        return false;
    }

    handler();

    return true;
}
