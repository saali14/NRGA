#version 440

uniform sampler3D TexUnit;
in vec3 MCposition;
 
out vec4 FragColor;

void main( void )
{
	vec3 stp = ( MCposition + 1.0 ) / 2.0; 
	FragColor = vec4( 0.8, 0.3, 0.1, 1.0 );
} 