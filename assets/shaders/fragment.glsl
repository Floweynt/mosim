#version 330 core

uniform vec3 diffuse_pos;  
uniform vec3 diffuse_color;
uniform vec3 ambient_color;
uniform float ambient_strength;

out vec4 FragColor;  
in vec4 out_color;
in vec3 out_normal;
in vec3 out_pos;
 

void main()
{
    float diffuse_strength = max(dot(out_normal, normalize(diffuse_pos - out_pos)), 0.0);
    FragColor = vec4((ambient_strength * ambient_color + diffuse_strength * diffuse_color) * vec3(out_color), 1);
    //FragColor = vec4(vec3(gl_FragCoord.z), 1.0);

}
