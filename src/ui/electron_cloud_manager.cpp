#include "ui/electron_cloud_manager.h"

void electron_cloud_manager::render_current_iso()
{
    if (result_to_render == nullptr)
    {
        return;
    }

    if (!dirty_status[mo_index])
    {
        return;
    }

    dirty_status[mo_index] = false;
    buffer_list[mo_index].isolevel_hf(*result_to_render, mo_index, CUBE_CORNER_1, CUBE_CORNER_2, SUBDIVISION_COUNT, ISOLEVEL);
}

void electron_cloud_manager::set_result(const hartree_fock_result* result)
{
    if (result == nullptr)
    {
        result_to_render = result;
        return;
    }

    result_to_render = result;
    size_t mo_count = result->coefficients.rows();
    dirty_status = std::vector<bool>(mo_count, true); // mark as dirty to trigger re-render

    effective_buffer_size = mo_count;
    if (effective_buffer_size > buffer_list.size())
    {
        buffer_list.resize(effective_buffer_size);
    }

    mo_index = 0; 
}

void electron_cloud_manager::render(gl::render_manager& render)
{
    if (result_to_render == nullptr)
    {
        return;
    }

    if (effective_buffer_size == 0)
    {
        return;
    }

    render_current_iso();
    buffer_list[mo_index].render(render, LIGHT_POS);
}

