#include "shaders.h"

#include <cstdio>
#include <fstream>
#include <vector>
#include <algorithm>


shaders::shaders()
{
	this->vertShaderId = 0;
	this->fragShaderId = 0;
	this->positionAttribute = 0;
	this->colourAttribute = 0;
	this->positionUniform = 0;
	this->accelerationAttribute = 0;
	this->velocityAttribute = 0;
	this->massAttribute = 0;
}

char* shaders::convertToString(const char *filename) 
{
	FILE *filePointer;
	char *content = NULL;

	int count = 0;

	if (filename != NULL) {
		filePointer = fopen(filename, "rt");

		if (filePointer != NULL) {

			fseek(filePointer, 0, SEEK_END);
			count = ftell(filePointer);
			rewind(filePointer);

			if (count > 0) {
				content = (char *)malloc(sizeof(char) * (count + 1));
				count = fread(content, sizeof(char), count, filePointer);
				content[count] = '\0';
			}
			fclose(filePointer);
		}
	}
	return content;
}

std::string shaders::readShadsersFile(const char* shaderFileName)
{
	std::string shaderCodes = this->convertToString(shaderFileName);
	
	return shaderCodes;
}

void shaders::compileShaders(const char* vertShaderFile, const char* fragShaderFile)
{
	GLenum      err;
	if ((err = glewInit()) != GLEW_OK)
		std::cout << glewGetErrorString(err) << std::endl;

	 this->vertShaderId = glCreateShader(GL_VERTEX_SHADER);
	 this->fragShaderId = glCreateShader(GL_FRAGMENT_SHADER);

	 // Read shaders
	 const char *vertShaderSrc = this->convertToString(vertShaderFile);
	 const char *fragShaderSrc = this->convertToString(fragShaderFile);

	 GLint result = GL_FALSE;
	 int logLength;

	 // Compile vertex shader
	 std::cout << "Compiling vertex shader." << std::endl;
	 glShaderSource(this->vertShaderId, 1, &vertShaderSrc, NULL);
	 glCompileShader(this->vertShaderId);

	 // Check vertex shader
	 glGetShaderiv(this->vertShaderId, GL_COMPILE_STATUS, &result);
	 glGetShaderiv(this->vertShaderId, GL_INFO_LOG_LENGTH, &logLength);
	 if (logLength > 0){
		 std::vector<char> VertexShaderErrorMessage(logLength + 1);
		 glGetShaderInfoLog(this->vertShaderId, logLength, NULL, &VertexShaderErrorMessage[0]);
		 printf("%s\n", &VertexShaderErrorMessage[0]);
	 }

	 // Compile fragment shader
	 std::cout << "Compiling fragment shader." << std::endl;
	 glShaderSource(this->fragShaderId, 1, &fragShaderSrc, NULL);
	 glCompileShader(this->fragShaderId);

	 // Check fragment shader
	 glGetShaderiv(this->fragShaderId, GL_COMPILE_STATUS, &result);
	 glGetShaderiv(this->fragShaderId, GL_INFO_LOG_LENGTH, &logLength);
	 if (logLength > 0){
		 std::vector<char> FragmentShaderErrorMessage(logLength + 1);
		 glGetShaderInfoLog(this->fragShaderId, logLength, NULL, &FragmentShaderErrorMessage[0]);
		 printf("%s\n", &FragmentShaderErrorMessage[0]);
	 }

	 std::cout << "Linking program" << std::endl;
	 this->shaderProgId = glCreateProgram();
	 glAttachShader(this->shaderProgId, this->vertShaderId);
	 glAttachShader(this->shaderProgId, this->fragShaderId);
	 glLinkProgram(this->shaderProgId);

	 glGetProgramiv(this->shaderProgId, GL_LINK_STATUS, &result);
	 glGetProgramiv(this->shaderProgId, GL_INFO_LOG_LENGTH, &logLength);
	 std::vector<char> programError((logLength > 1) ? logLength : 1);
	 glGetProgramInfoLog(this->shaderProgId, logLength, NULL, &programError[0]);
	 std::cout << &programError[0] << std::endl;


	 // Detaching and deleting shaders from program 
	 glDetachShader(this->shaderProgId, this->vertShaderId);
	 glDetachShader(this->shaderProgId, this->fragShaderId);

	 glDeleteShader(this->vertShaderId);
	 glDeleteShader(this->fragShaderId);
}


GLuint shaders::getProgShaderId()
{
	return this->shaderProgId;
}

GLuint shaders::getPositionAttribute()
{
	return this->positionAttribute;
}

GLuint shaders::getColorAttribute()
{
	return this->colourAttribute;
}

GLuint shaders::getNormalAttribute()
{
	return this->normalAttribute;
}

GLuint shaders::getAccelerationAttribute()
{
	return this->accelerationAttribute;
}

GLuint shaders::getVelocityAttribute()
{
	return this->velocityAttribute;
}

GLuint shaders::getMassAttribute()
{
	return this->massAttribute;
}

GLuint shaders::getPositionUniform()
{
	return this->positionUniform;
}

GLuint shaders::getProjectionMatrixUniform()
{
	return this->MPuniform;
}

GLuint shaders::getModelViewMatrixUniform()
{
	return this->MVuniform;
}


void shaders::setPositionAttribute(GLuint positionAttrib)
{
	this->positionAttribute = positionAttrib;
}

void shaders::setColorAttribute(GLuint colorAttrib)
{
	this->colourAttribute = colorAttrib;
}

void shaders::setNormalAttribute(GLuint normalAttrib)
{
	this->normalAttribute = normalAttrib;
}

void shaders::setAccelerationAttribute(GLuint accelerationAttrib)
{
	this->accelerationAttribute = accelerationAttrib;
}

void shaders::setVelocityAttribute(GLuint velocityAttrib)
{
	this->velocityAttribute = velocityAttrib;
}

void shaders::setMassAttribute(GLuint massAttrib)
{
	this->massAttribute = massAttrib;
}

void shaders::setPositionUniform(GLuint posUniform)
{
	 this->positionUniform = positionUniform;
}

void shaders::setModelViewMatrixUniform(GLuint mv)
{
	this->MVuniform = mv;
}

void shaders::setProjectionMatrixUniform(GLuint mp)
{
	this->MPuniform = mp;
}

shaders::~shaders()
{
}
