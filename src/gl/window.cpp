// cSpell:ignore glfw glew
#include <GL/glew.h>
// dont reorder
#include "gl/window.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <fmt/core.h>
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/glm.hpp>
#include <iostream>
#include <stdexcept>

gl::window::window(int width, int height, const char* name)
{
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);

    GLFWwindow* win = glfwCreateWindow(width, height, name, nullptr, nullptr);

    if (!win)
    {
        throw std::runtime_error("Failed to open GLFW window");
    }

    glfwMakeContextCurrent(win);
    auto glew_ret = glewInit();
    if (glew_ret != GLEW_OK)
    {
        throw std::runtime_error(fmt::format("Failed to initialize GLEW: {}", (const char*)glewGetErrorString(glew_ret)));
    }

    glfwMakeContextCurrent(win);
    glfwSetWindowUserPointer(win, this);

    // setup callbacks for the event handler functions
    glfwSetKeyCallback(
        win, +[](GLFWwindow* win, int key, int scancode, int action, int mods) {
            auto* that = (window*)glfwGetWindowUserPointer(win);
            if (action == GLFW_PRESS)
            {
                that->key_pressed(key, scancode, mods);
            }
            else if (action == GLFW_RELEASE)
            {
                that->key_released(key, scancode, mods);
            }
        });

    glfwSetMouseButtonCallback(
        win, +[](GLFWwindow* win, int button, int action, int mods) {
            auto* that = (window*)glfwGetWindowUserPointer(win);
            if (action == GLFW_PRESS)
            {
                that->mouse_button_pressed(button, mods);
            }
            else if (action == GLFW_RELEASE)
            {
                that->mouse_button_released(button, mods);
            }
        });

    glfwSetScrollCallback(
        win, +[](GLFWwindow* win, double x, double y) {
            auto* that = (window*)glfwGetWindowUserPointer(win);
            that->mouse_scrolled(x, y);
        });

    glfwSetCursorPosCallback(
        win, +[](GLFWwindow* win, double mouse_x, double mouse_y) {
            auto* that = (window*)glfwGetWindowUserPointer(win);
            that->mouse_moved(mouse_x, mouse_y);
        });

    glfwSwapInterval(1);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    this->win = win;
}

gl::window::~window()
{
    glfwDestroyWindow(win);
}

void gl::window::close() { glfwSetWindowShouldClose(win, GLFW_TRUE); }

void gl::window::run()
{
    init();

    while (!glfwWindowShouldClose(win))
    {
        GLint window_width = 0;
        GLint window_height = 0;
        glfwGetWindowSize(win, &window_width, &window_height);
        glViewport(0, 0, window_width, window_height);

        glClearColor(bg[0], bg[1], bg[2], bg[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        render(glm::perspective(fov, (double)window_width / (double)window_height, z_near, z_far), {window_width, window_height});

        glfwSwapBuffers(win);
        glfwPollEvents();
    }
}
