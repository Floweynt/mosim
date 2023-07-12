#pragma once

#include <glm/glm.hpp>
#include <string_view>

#define TEXT_LEFT 0
#define TEXT_TOP 0
#define TEXT_CENTER 1
#define TEXT_RIGHT 2
#define TEXT_BOTTOM 2

namespace gl
{
    void render_text(const std::string_view& str, glm::vec2 pos, float scale = 1, glm::vec4 color = {1, 1, 1, 1}, int horiz = TEXT_LEFT,
                     int vert = TEXT_TOP);
} // namespace gl
