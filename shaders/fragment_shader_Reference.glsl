#version 440

in vec4 FrontColor;
in vec4 BackColor;

out vec4 FragColor;


void main(void)
{
	FragColor = BackColor;
}