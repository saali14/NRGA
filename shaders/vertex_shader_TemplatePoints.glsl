#version 440

in vec3  position;
in vec3  colour;
in vec3  acceleration;
in vec3  velocity;
in float mass;


struct LightInfo {
	vec4 Position;	// Light position in eye coords.
	vec3 La;		// Ambient light intensity
	vec3 Ld;		// Diffuse light intensity
	vec3 Ls;		// Specular light intensity
};

uniform LightInfo		Light;
uniform mat4			MV;
uniform mat4			MP;
uniform sampler1D		ColorMap;

out vec3 MCposition; 
out vec3 clr;

void main(void)
{
	clr = colour/255.0;
	MCposition  = position;
	gl_Position = MP * MV * vec4(position, 1.0);
} 


