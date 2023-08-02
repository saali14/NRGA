#version 440

in vec3 FrontColor;
in vec3 BackColor;

out vec4 FragColor;


void main(void)
{
	FragColor = vec4(BackColor, 1.0);
}