#version 440

uniform sampler3D TexUnit;
in vec3 MCposition; 
in vec3 clr;

out vec4 FragColor;

void main(void)
{
	vec3 stp = ( MCposition + 1.0 ) / 2.0;
	FragColor = vec4( clr , 1.0 );
} 
