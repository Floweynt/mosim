#ifndef __WINDOW_H__
#define __WINDOW_H__

#include <GLFW/glfw3.h>
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/ext/vector_float4.hpp>
#include <iostream>
#include <stdexcept>

namespace gl
{
    using namespace glm;
    class window
    {
        GLFWwindow* win;

    protected:
        double fov;
        double z_near;
        double z_far;

        glm::vec4 bg;
        virtual void init() = 0;
        virtual void render(const mat4& proj, glm::uvec2 screen_size) = 0;
        virtual void key_pressed(int key, int scancode, int mods) = 0;
        virtual void key_released(int key, int scancode, int mods) = 0;
        virtual void mouse_moved(double xpos, double ypos) = 0;
        virtual void mouse_button_pressed(int button, int mods) = 0;
        virtual void mouse_button_released(int button, int mods) = 0;
        virtual void mouse_scrolled(double x_off, double y_off) = 0;

        void close();

    public:
        window(int width, int height, const char* name);
        virtual ~window();
        void run();
    };
} // namespace gl

#endif
