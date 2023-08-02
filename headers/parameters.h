/*
 * GA_settings.hpp
 *
 *  Created on: 16.09.2016
 *      Author: golyanik, ali
 * Modified on: 28.07.2017
 * 
 */

#ifndef PARAMETER_HEADERS
#define PARAMETER_HEADERS

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

#ifdef _WIN32
#include <libconfig.h>
#endif

#ifdef __linux__
#include <libconfig.h++>
#endif

using std::string;
using namespace libconfig;

template<typename T>
struct parameters
{
	int i;
	int p_Reps;
	
	unsigned int color[3];

	string  in_template,
		in_reference,
		in_template_subsample,
		in_reference_subsample,
		in_Transformation,
		out_registered;
			
	int 	p_numPtsTemplate,
		p_numPtsReference;
			
				 
	T       p_alpha_normalized,
		p_delta,
		p_convergeErrorRario;
			
	int 	p_max_iterations, 
		p_useHuber;
    			
	T       p_window_width,
		p_window_height,
		p_window_depth,
		p_left_bound,
		p_right_bound,
		p_top_bound,
		p_bottom_bound;

	T       p_dimension,
		p_theta_X,
		p_theta_Y,
		p_theta_Z,
		p_translation_x,
		p_translation_y,
		p_translation_z,
		p_scale_xyz;

	int     p_color_output,
		p_addnoise_flag,
		p_normalize_flag,
		p_mass_normalize_flag;
			
	std::string  p_noise_types;
			
	T       p_noiseG_variance,
		p_noiseG_mean,
		p_noiseU_Xmin,
		p_noiseU_Xmax,
		p_noiseU_Ymin,
		p_noiseU_Ymax,
		p_noiseU_Zmin,
		p_noiseU_Zmax,
		p_noise_percentage,
		p_tol;
		    
	// Gravitational PARAMETERS
	T       p_forceX_Init;
	T       p_forceY_Init;
	T       p_forceZ_Init;
	T       p_accX_Init;
	T       p_accY_Init;
	T       p_accZ_Init;
	T       p_mass;
	T       p_massT;
	T       p_massR;
	T       p_G;
	T       p_GPE;
	T       p_velX_Init;
	T       p_velY_Init;
	T       p_velZ_Init;
	T       p_epsilon;
	T       p_ballRadius;
	T       p_dilatation;
	T       p_spaceCurvature;
	T       p_timeStep;

	// MASS SPRINGS PARAMETERS
	T       p_spring_Constraint;
	T       p_damping_Constraint;
	T       p_stretch_Constraint;
	T       p_bend_Constraint;
	T       p_spring_Mass;
	T       p_rotation_max;

	// NEAREST NEIGHBOURS PARAMETERS
	int  k_T;
	int  k_R;
	int  dimT;
	int  dimR;

	// VISUALIZATION PARAMETERS
	bool    useMaterialColor;
	bool    useShadersT;
	bool    useShadersR;
	bool    useShadersRMesh;
	bool    useShadersRPoints;
	bool    useShadersTMesh;
	bool    useShadersTPoints;
	
	float   bgColor[4];
        float   transparency = 0.8;
	float   LightPositionT[4];
	float   LDeffuseColorT[3];
	float   LSpecularColorT[3];
	float   LAmbientColorT[3];
	float   MAmbientT[3];
	float   MDiffuseT[3];
	float   MSpecularT[3];
	float   MShinynessT;
	float   LightPositionR[4];
	float   LDeffuseColorR[3];
	float   LSpecularColorR[3];
	float   LAmbientColorR[3];
	float   MAmbientR[3];
	float   MDiffuseR[3];
	float   MSpecularR[3];
	float   MShinynessR;
	float   pointSize;	

	string  shaderPath_TMeshVert;
	string  shaderPath_RMeshVert;
	string  shaderPath_TMeshFrag;
	string  shaderPath_RMeshFrag;
	string  shaderPath_TPointsVert;
	string  shaderPath_RPointsVert;
	string  shaderPath_TPointsFrag;
	string  shaderPath_RPointsFrag;

	// EVALUATIONS PARAMETERS
	void load_settings(string config_file);
	void dump_settings();
};


template<class T> void parameters<T>::load_settings(string config_file){

	libconfig::Config cfg;

	try
	{
	  cfg.readFile(config_file.c_str());
	  libconfig::Setting &root = cfg.getRoot();

	  int color_flag                        = root["color_flag"];
	  libconfig::Setting &default_color 	= root["p_color_default"];
	  libconfig::Setting &template_color 	= root["p_color_template"];
	  libconfig::Setting &refernce_color 	= root["p_color_reference"];	
	  p_color_output                        = root["p_color_output"];
	
	  /*
	  if(p_color_output == 0)
	  {
	    for (uint i = 0; i < 3; i++)
	      color[i] = rgb[i];
	  }
	  else if(p_color_output == 1)
	  {
	    for (uint i = 0; i < 3; i++)
	      color[i] = rgb[i];
	  }
	  else if(p_color_output == 2)
	  {
	    for (uint i = 0; i < 3; i++)
	      color[i] = rgb[i];
	  }
	  else
	  {
	    for (uint i = 0; i < 3; i++)
	      color[i] = rgb[i];
	  }
	  */
	  
	  string str_param1 		= root["in_template"];
	  string str_param2 		= root["in_reference"];
	  string str_param3 		= root["in_template_subsample"];
	  string str_param4 		= root["in_reference_subsample"];
	  string str_param5         	= root["in_Transformation"];
	  string str_param6 		= root["out_registered"];
	  
	  string str_param7 		= root["in_shaderPath_TMeshVert"];
	  string str_param8 		= root["in_shaderPath_RMeshVert"];
	  string str_param9 		= root["in_shaderPath_TMeshFrag"];
	  string str_param10 		= root["in_shaderPath_RMeshFrag"];
	  string str_param11 		= root["in_shaderPath_TPointsVert"];
	  string str_param12 		= root["in_shaderPath_RPointsVert"];
	  string str_param13 		= root["in_shaderPath_TPointsFrag"];
	  string str_param14 		= root["in_shaderPath_RPointsFrag"];
	
	  in_template               	= str_param1;
	  in_reference              	= str_param2;
	  in_template_subsample		= str_param3;
	  in_reference_subsample	= str_param4;
	  in_Transformation         	= str_param5;
	  out_registered            	= str_param6;
	
	  shaderPath_TMeshVert		= str_param7;
	  shaderPath_RMeshVert		= str_param8;
	  shaderPath_TMeshFrag		= str_param9;
	  shaderPath_RMeshFrag		= str_param10;
	  shaderPath_TPointsVert	= str_param11;
	  shaderPath_RPointsVert	= str_param12;
	  shaderPath_TPointsFrag	= str_param13;
	  shaderPath_RPointsFrag	= str_param14;
	
	  p_window_width            	= root["p_window_width"];
	  p_window_height           	= root["p_window_height"];
	  p_window_depth            	= root["p_window_depth"];
	  p_left_bound              	= root["p_left_bound"];
	  p_right_bound             	= root["p_right_bound"];
	  p_top_bound               	= root["p_top_bound"];
	  p_bottom_bound            	= root["p_bottom_bound"];

	  p_addnoise_flag           	= root["p_addnoise_flag"];
	  p_normalize_flag          	= root["p_normalize_flag"];
	  p_mass_normalize_flag		= root["p_mass_normalize_flag"];
	  
	  string str_noisetype      	= root["p_noise_type"];
	  p_noise_types             	= str_noisetype;
	  p_noiseG_variance         	= root["p_noiseG_sigma"];
	  p_noiseG_mean             	= root["p_noiseG_mean"];
	  p_noiseU_Xmin             	= root["p_noiseU_Xmin"];
	  p_noiseU_Xmax             	= root["p_noiseU_Xmin"];
	  p_noiseU_Ymin             	= root["p_noiseU_Ymin"];
	  p_noiseU_Ymax             	= root["p_noiseU_Ymax"];
	  p_noiseU_Zmin             	= root["p_noiseU_Zmin"];
	  p_noiseU_Zmax             	= root["p_noiseU_Zmax"];
	  p_noise_percentage		= root["p_noise_percentage"];
	  p_tol                     	= root["p_tol"];
	  p_Reps			= root["p_Reps"];
	  p_max_iterations		= root["p_max_iterations"];
	  
	  // Gravitational PARAMETERS
	  p_dimension               	= root["p_dimension"];
	  p_G                       	= root["p_G"];
	  p_mass                    	= root["p_mass"];
	  p_massT                   	= root["p_MT"];
	  p_massR                   	= root["p_MR"];
	  p_epsilon                 	= root["p_epsilon"];
	  p_GPE                     	= root["p_GPE"];
	  p_ballRadius              	= root["p_ballRadius"];
	  p_dilatation              	= root["p_dilatation"];
	  p_spaceCurvature          	= root["p_spaceCurvature"];
	  p_velX_Init               	= root["p_vel_X"];
	  p_velY_Init               	= root["p_vel_Y"];
	  p_velZ_Init               	= root["p_vel_Z"];
	  p_timeStep                	= root["p_timeStep"];
	  p_theta_X			= root["p_theta_X"];
	  p_theta_Y			= root["p_theta_Y"];
	  p_theta_Z			= root["p_theta_Z"];
	  p_translation_x           	= root["p_translation_x"];
	  p_translation_y           	= root["p_translation_y"];
	  p_translation_z           	= root["p_translation_z"];
	  p_scale_xyz               	= root["p_scale_xyz"];

	  // MASS SPRINGS PARAMETERS
	  p_spring_Constraint		= root["p_spring_Constraint"];
	  p_damping_Constraint		= root["p_damping_Constraint"];
	  p_stretch_Constraint		= root["p_stretch_Constraint"];
	  p_bend_Constraint         	= root["p_bend_Constraint"];
	  p_spring_Mass             	= root["p_spring_Mass"];
	  p_rotation_max            	= root["p_rotation_max"];

	  // NEAREST NEIGHBOURS PARAMETERS
	  k_T                       	= root["k_T"];			;
	  k_R                       	= root["k_R"];
	  dimT                      	= root["dimT"];
	  dimR                      	= root["dimR"];

	  // VISUALIZATION PARAMETERS
	  useMaterialColor          	= root["useMaterialColor"];
	  useShadersT               	= root["useShadersT"];
	  useShadersR               	= root["useShadersR"];

	  libconfig::Setting &bgc   	= root["bgColor"];
	  
	  libconfig::Setting &LposT  	= root["LightPositionT"];
	  libconfig::Setting &LDclrT 	= root["LDeffuseColorT"];
	  libconfig::Setting &LSclrT 	= root["LSpecularColorT"];
	  libconfig::Setting &LAclrT 	= root["LAmbientColorT"];
	  libconfig::Setting &MDclrT 	= root["MDiffuseT"];
	  libconfig::Setting &MSclrT 	= root["MSpecularT"];
	  libconfig::Setting &MAclrT 	= root["MAmbientT"];
	  
	  libconfig::Setting &LposR  	= root["LightPositionR"];
	  libconfig::Setting &LDclrR 	= root["LDeffuseColorR"];
	  libconfig::Setting &LSclrR 	= root["LSpecularColorR"];
	  libconfig::Setting &LAclrR 	= root["LAmbientColorR"];
	  libconfig::Setting &MDclrR 	= root["MDiffuseR"];
	  libconfig::Setting &MSclrR 	= root["MSpecularR"];
	  libconfig::Setting &MAclrR 	= root["MAmbientR"];

	  for (uint i = 0; i < 4; i++)
	    bgColor[i] = bgc[i];

	  for (uint i = 0; i < 4; i++)
	  {
	    LightPositionT[i] = LposT[i];
	    LightPositionR[i] = LposR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    LDeffuseColorT[i] = LDclrT[i];
	    LDeffuseColorR[i] = LDclrR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    LSpecularColorT[i] = LSclrT[i];
	    LSpecularColorR[i] = LSclrR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    LAmbientColorT[i] = LAclrT[i];
	    LAmbientColorR[i] = LAclrR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    MDiffuseT[i] = MDclrT[i];
	    MDiffuseR[i] = MDclrR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    MSpecularT[i] = MSclrT[i];
	    MSpecularR[i] = MSclrR[i];
	  }

	  for (uint i = 0; i < 3; i++)
	  {
	    MAmbientT[i] = MAclrT[i];
	    MAmbientR[i] = MAclrR[i];
	  }

	  MShinynessT                = root["MShinynessT"];
	  MShinynessR                = root["MShinynessR"];
	  pointSize                 = root["pointSize"];
	  useShadersRMesh           = root["useShadersRMesh"];
	  useShadersRPoints         = root["useShadersRPoints"];
	  useShadersTMesh           = root["useShadersTMesh"];
	  useShadersTPoints         = root["useShadersTPoints"];
	}
	catch(const libconfig::FileIOException &fioex){
		std::cerr << "-> config file not found" << std::endl;
		exit(EXIT_FAILURE);
	}
	catch( const libconfig::SettingTypeException &ste ){
	  std::cerr << "-> failed to read settings file" << std::endl;
	  std::cerr << "  -> possible reason: type mismatch (e.g. float parameters should contain a \".\")" << std::endl;
	  exit(EXIT_FAILURE);  
	}
	catch (const libconfig::SettingNotFoundException &snf){
	  std::cerr << "-> failed to read settings file" << std::endl;
	  std::cerr << "  -> a setting accidently deleted from the settings file" << std::endl;
	  exit(EXIT_FAILURE);
	}
	catch(...){
	  std::cerr << "-> failed to read settings file" << std::endl;
	  std::cerr << "  -> possible reason: the file corrupted" << std::endl;
	  exit(EXIT_FAILURE);
	}
}

//dump settings formatted
template<class T>  void parameters<T>::dump_settings(){
}


#endif
