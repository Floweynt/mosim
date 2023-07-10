#pragma once

#include <glm/glm.hpp>
#include <string_view>

namespace gl
{
    void draw_text(const std::string_view& str, glm::vec2 pos, float scale = 1, glm::vec4 color = {1, 1, 1, 1}, glm::vec4 bg = {0, 0, 0 ,0});
} // namespace gl
