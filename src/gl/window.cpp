// cSpell:ignore glfw glew
#include <GL/glew.h>
// dont reorder
#include "gl/window.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/glm.hpp>
#include <iostream>

gl::window::window(int width, int height, const char* name)
{
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

    GLFWwindow* win = glfwCreateWindow(width, height, name, nullptr, nullptr);

    if (!win)
        throw std::runtime_error("Failed to open GLFW window");

    glfwMakeContextCurrent(win);
    auto glew_ret = glewInit();
    if (glew_ret != GLEW_OK)
    {
        std::cerr << "Failed to initialize GLEW: msg=" << glewGetErrorString(glew_ret) << "\n";
        std::exit(-1);
    }

    glfwMakeContextCurrent(win);
    glfwSetWindowUserPointer(win, this);
    glfwSetKeyCallback(
        win, +[](GLFWwindow* w, int key, int scancode, int action, int mods) {
            window* that = (window*)glfwGetWindowUserPointer(w);
            that->controls(key, scancode, action, mods);
        });
    glfwSetMouseButtonCallback(
        win, +[](GLFWwindow* w, int button, int action, int mods) {
            window* that = (window*)glfwGetWindowUserPointer(w);
            that->mause_button(button, action, mods);
        });
    glfwSetScrollCallback(
        win, +[](GLFWwindow* w, double x, double y) {
            window* that = (window*)glfwGetWindowUserPointer(w);
            that->mause_scroll(x, y);
        });
    glfwSetCursorPosCallback(
        win, +[](GLFWwindow* w, double mouse_x, double mouse_y) {
            window* that = (window*)glfwGetWindowUserPointer(w);
            that->mause(mouse_x, mouse_y);
        });
    glfwSwapInterval(1);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    this->win = win;
}

gl::window::~window() { glfwDestroyWindow(win); }

void gl::window::close() { glfwSetWindowShouldClose(win, GL_TRUE); }

void gl::window::run()
{
    init();

    while (!glfwWindowShouldClose(win))
    {
        GLint window_width, window_height;
        glfwGetWindowSize(win, &window_width, &window_height);
        glViewport(0, 0, window_width, window_height);

        glClearColor(bg[0], bg[1], bg[2], bg[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        render(glm::perspective(fov, (double)window_width / (double)window_height, z_near, z_far));

        glfwSwapBuffers(win);
        glfwPollEvents();
    }
}
