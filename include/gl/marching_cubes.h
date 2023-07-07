#pragma once
#include "gl/glstate.h"
#include "hf/hf.h"
#include "gl/render.h"

class hf_isosurface
{
    gl::vertex_array_object vao;
    gl::buffer_object vbo;
    size_t nverts;

public:
    
    void draw(gl::render_manager& render);
    void isolevel_hf(const hartree_fock_result& result, uint32_t which_mo, glm::vec3 start, glm::vec3 end, size_t detail);
};

