#version 440

in vec3 position;
in vec3 colour;

struct LightInfo {
	vec4 Position;	// Light position in eye coords.
	vec3 La;		// Ambient light intensity
	vec3 Ld;		// Diffuse light intensity
	vec3 Ls;		// Specular light intensity
};


uniform LightInfo Light;
uniform mat4 MV;
uniform mat4 MP;

out vec3 MCposition; 

void main( void )
{
	MCposition  = position;
	gl_Position = MP * MV * vec4(position, 1.0);
} 