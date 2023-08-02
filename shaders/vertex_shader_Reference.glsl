#version 440

in vec3 position;
in vec3 colour;
in vec3 normal;

out vec4 FrontColor;
out vec4 BackColor;

struct LightInfo {
	vec4 Position;		// Light position in eye coords.
	vec3 La;		// Ambient light intensity
	vec3 Ld;		// Diffuse light intensity
	vec3 Ls;		// Specular light intensity
};

struct MaterialInfo {
	vec3 Ka;		// Ambient reflectivity
	vec3 Kd;		// Diffuse reflectivity
	vec3 Ks;		// Specular reflectivity
	float Shininess; 	// Specular shininess factor
};

uniform LightInfo Light;
uniform MaterialInfo Material;
uniform float p;
uniform mat4 MV;
uniform mat4 MP;

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
		
        
        float distanceToLight = length(normalize(Light.Position - epos));
        float attenuation = 1.0 / (1.0 + 0.1 * pow(distanceToLight, 2));

	vec3 linearColor = ambient +  (attenuation * (diffuse + spec));
	
	vec3 gamma = vec3(1.0/2.2);
        vec3 finalColor = vec3(pow(linearColor, gamma));
        
        return linearColor;
}

void main(void)
{
	vec3 tnorm = normalize(mat3(transpose(inverse(MV))) * normal);
	vec4 eyeCoords = MV * vec4(position, 1.0);
	FrontColor = vec4(phongModel(eyeCoords, tnorm), 0.8);
	BackColor = vec4(phongModel(eyeCoords, -tnorm), 0.8);
	gl_Position = MP * MV * vec4(position, 1.0);
}
