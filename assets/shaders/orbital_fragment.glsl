#version 330 core

uniform vec3 diffuse_pos;  
uniform vec3 diffuse_color;
uniform vec3 ambient_color;
uniform float ambient_strength;

out vec4 FragColor;  

in vec4 geom_color;
in vec3 geom_normal;
in vec3 geom_pos;
 
void main()
{
    float diffuse_strength = max(dot(geom_normal, normalize(diffuse_pos - geom_pos)), 0.0);
    FragColor = vec4((ambient_strength * ambient_color + diffuse_strength * diffuse_color) * vec3(geom_color), 1);
}
