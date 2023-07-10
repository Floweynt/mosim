#include <GL/glew.h>
//
#include <GL/gl.h>
#include <stdexcept>
//
#include "gl/font.h"
#include "singleton.h"
#define GLT_IMPLEMENTATION
#include "gl/gltext.h"

namespace gl
{
    class glfont_handle : public singleton<glfont_handle>
    {
    protected:
        glfont_handle()
        {
            if (!gltInit())
            {
                throw std::runtime_error("Failed to initialize glText\n");
            }
        }

        ~glfont_handle()
        {
            gltTerminate();
        }
    public:
    };

    void draw_text(const std::string_view& str, glm::vec2 pos, float scale, glm::vec4 color, glm::vec4 bg)
    {
        glfont_handle::get_instance();

        GLTtext* text = gltCreateText();
        gltSetText(text, str.data());

        gltBeginDraw();
        gltColor(color.r, color.g, color.b, color.a);
        gltDrawText2D(text, pos.x, pos.y, scale);
        gltEndDraw();

        gltDeleteText(text);
    }

} // namespace gl
