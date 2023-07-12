#pragma once

#include <glm/glm.hpp>

inline static constexpr auto KB_UI_SIZE = glm::vec2(500, 50);
inline static constexpr auto CONFIG_UI_SIZE = KB_UI_SIZE;
inline static constexpr auto LABEL_OFFSET = 10;
inline static constexpr auto MAIN_KB_OFFSET = 300;
inline static constexpr auto ALT_KB_OFFSET = 400;
inline static constexpr auto KB_DISPLAY_SIZE = glm::vec2(90, 30);
inline static constexpr glm::dvec3 CUBE_CORNER_1 = {-6, -6, -6};
inline static constexpr glm::dvec3 CUBE_CORNER_2 = {6, 6, 6};
inline static constexpr auto SUBDIVISION_COUNT = 20;
inline static constexpr auto ISOLEVEL = 0.05;
inline static constexpr glm::vec3 LIGHT_POS = {10, 10, 10};
inline static constexpr auto SPHERE_SUBDIV_A = 0.02;
inline static constexpr auto SPHERE_SUBDIV_B = 0.04;
inline static constexpr uint32_t BOX_COLOR = 0xc8c8ff;
inline static constexpr auto ATOM_SPHERE_RADIUS = 0.1;
inline static constexpr auto MOUSE_MOVE_FACTOR = 0.005;
inline static constexpr auto MOUSE_ZOOM_FACTOR = 0.693;
inline static constexpr int64_t MARGIN_BOTTOM = 20;

inline static constexpr glm::vec3 PASSIVE_COLOR = {0.5, 0.5, 0.5};
inline static constexpr glm::vec3 ACTIVE_COLOR = {0.7, 0.7, 0.7};
inline static constexpr glm::vec3 CONFIG_BG_COLOR = {0.4, 0.4, 0.4};
inline static constexpr glm::vec3 TEXT_COLOR = {1, 1, 1};

