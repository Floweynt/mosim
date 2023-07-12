#include <GL/glew.h>
//
#include <stdexcept>
//
#include "gl/font.h"
#include "singleton.h"
#define GLT_IMPLEMENTATION
#include "gl/gltext.h"

namespace gl
{
    class gl_text_handle : public singleton<gl_text_handle>
    {
    protected:
        gl_text_handle()
        {
            if (!gltInit())
            {
                throw std::runtime_error("Failed to initialize glText\n");
            }
        }

        ~gl_text_handle() { gltTerminate(); }

    public:
    };

    void render_text(const std::string_view& str, glm::vec2 pos, float scale, glm::vec4 color, int horiz, int vert)
    {
        gl_text_handle::get_instance();

        GLTtext* text = gltCreateText();
        gltSetText(text, str.data());

        gltBeginDraw();
        gltColor(color.r, color.g, color.b, color.a);

        gltDrawText2DAligned(text, pos.x, pos.y, scale, horiz, vert);

        gltEndDraw();

        gltDeleteText(text);
    }

} // namespace gl
