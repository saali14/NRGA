#version 440

in vec3  position;
in vec3  colour;
in vec3  acceleration;
in vec3  velocity;
in float mass;
in vec3  normal;

out vec3 FrontColor;
out vec3 BackColor;

struct LightInfo {
	vec4 Position;	// Light position in eye coords.
	vec3 La;		// Ambient light intensity
	vec3 Ld;		// Diffuse light intensity
	vec3 Ls;		// Specular light intensity
};

struct MaterialInfo {
	vec3 Ka;		// Ambient reflectivity
	vec3 Kd;		// Diffuse reflectivity
	vec3 Ks;		// Specular reflectivity
	float Shininess; // Specular shininess factor
};

uniform LightInfo		Light;
uniform MaterialInfo		Material;
uniform vec2			p;
uniform mat4			MV;
uniform mat4			MP;
uniform sampler1D		ColorMap;

vec3 phongModel(vec4 epos, vec3 norm)
{
	vec3 s = normalize(vec3(Light.Position - epos));
	vec3 v = normalize(-epos.xyz);
	vec3 r = reflect(-s, norm);
	vec3 ambient = Light.La * Material.Ka;

	float sDotN = max(dot(s, norm), 0.0);
	vec3 diffuse = Light.Ld * Material.Kd * sDotN;
	
	vec3 spec = vec3(0.0);

	if (sDotN > 0.0)
		spec = Light.Ls * Material.Ks * pow(max(dot(r, v), 0.0), Material.Shininess);
	return ambient + diffuse + spec;
}

void main(void)
{
	vec3 tnorm = normalize(mat3(transpose(inverse(MV))) * normal);
	vec4 eyeCoords = MV * vec4(position, 1.0);
	FrontColor = phongModel(eyeCoords, tnorm);
	BackColor = phongModel(eyeCoords, -tnorm);
	gl_Position = MP * MV * vec4(position, 1.0);
}
