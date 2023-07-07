#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <iostream>

[[gnu::constructor]] void init()
{
    glfwSetErrorCallback(+[](int, const char* msg )
            { std::cerr << msg << '\n'; });

    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        std::exit(-1);
    }
}

[[gnu::destructor]] void exit() { glfwTerminate(); }
