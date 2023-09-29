#pragma once

#include "gl/orbital_render.h"
#include "hf/hf.h"
#include "ui_config.h"

class electron_cloud_manager
{
    const hartree_fock_result* result_to_render{};
    std::vector<bool> dirty_status;
    std::vector<hf_isosurface> buffer_list;
    size_t effective_buffer_size{};
    size_t mo_index{};
    double isolevel = ISOLEVEL;

    void render_current_iso();

public:
    constexpr void next_mo() { mo_index = (mo_index + 1) % effective_buffer_size; }
    constexpr void prev_mo() { mo_index = (mo_index + effective_buffer_size - 1) % effective_buffer_size; }
    [[nodiscard]] constexpr auto get_index() const { return mo_index; }
    constexpr void set_mo(size_t index)
    {
        if (index < effective_buffer_size)
        {
            mo_index = index;
        }
    }

    void set_result(const hartree_fock_result* result);
    void render(gl::render_manager& render);
    constexpr void regenerate() { dirty_status = std::vector<bool>(effective_buffer_size, true); }

    void add_isolevel(double value);
};

