#include "nfd.h"
#include <GL/glew.h>
//
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <iostream>

[[gnu::constructor]] void init()
{
    glfwSetErrorCallback(+[](int, const char* msg) { std::cerr << msg << '\n'; });

    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        std::exit(-1);
    }

    NFD_Init();
    if (NFD_Init() != NFD_OKAY)
    {
        throw std::runtime_error("failed to init NFD");
    }
}

[[gnu::destructor]] void exit()
{
    NFD_Quit();

    glfwTerminate();
}
