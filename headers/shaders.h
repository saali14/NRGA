#ifndef GLSL_SHADERS_NRGA
#define GLSL_SHADERS_NRGA

#pragma once
#include <iostream>
#include <string>

#define GL_STATIC
#include <GL/glew.h>

class shaders
{
public:
	shaders();
	~shaders();
	char* convertToString(const char *filename);
	std::string readShadsersFile(const char* shaderFileName);
	void compileShaders(const char* vertShaderFile, const char* fragShaderFile);

	GLuint getProgShaderId();
	GLuint getPositionAttribute();
	GLuint getColorAttribute();
	GLuint getNormalAttribute();
	GLuint getAccelerationAttribute();
	GLuint getVelocityAttribute();
	GLuint getMassAttribute();
	GLuint getPositionUniform();
	GLuint getProjectionMatrixUniform();
	GLuint getModelViewMatrixUniform();
	
	void setPositionAttribute(GLuint);
	void setColorAttribute(GLuint);
	void setNormalAttribute(GLuint);
	void setAccelerationAttribute(GLuint);
	void setVelocityAttribute(GLuint);
	void setMassAttribute(GLuint);
	void setPositionUniform(GLuint);
	void setProjectionMatrixUniform(GLuint);
	void setModelViewMatrixUniform(GLuint);

private:
	GLuint positionUniform, MPuniform, MVuniform; 
	GLuint colourAttribute, positionAttribute, normalAttribute, accelerationAttribute, velocityAttribute, massAttribute;
	GLuint vertShaderId, fragShaderId, shaderProgId;
};

#endif // GLSL_SHADERS_NRGA
