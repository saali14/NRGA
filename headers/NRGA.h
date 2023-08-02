#ifndef _NRGA_
#define _NRGA_

/*
 *  @Author : Sk Aziz Ali
 *  @Place  : DFKI GmbH, Kaiserslautern, Germany
 *  @Date   : 22nd October, 2016
*/

#pragma once

#include <headers/gnuplot-iostream.h>
#include "parameters.h"
#include "KD_Tree.h"

template<typename T>
class NRGA
{
public:
	NRGA();
	~NRGA();

	T meanDeviationTR;
	T RMSE_TR;

	Gnuplot gp;
	std::vector<T> GPE_values;

	parameters<T>* param;
	void		NRGA_setParameters(parameters<T>* param);
	parameters<T>* 	NRGA_getParameters();
	
	void NRGA_setT(shape_point<T>* Tpt);
	shape_point<T>* NRGA_getT();
	void NRGA_initShapePoints_T3D(shape_point<T> &);
	void NRGA_initShapePoints_R3D(shape_point<T> &);
	void NRGA_initShapePoints_T4D(shape_point<T> &);
	void NRGA_initShapePoints_R4D(shape_point<T> &);
	
	
	point3D<T>* NRGA_getShapePoints_T3DPrev();
	point4D<T>* NRGA_getShapePoints_T4DPrev();
	point3D<T>* NRGA_getShapePoints_T3D();
	point3D<T>* NRGA_getShapePoints_R3D();
	point4D<T>* NRGA_getShapePoints_T4D();
	point4D<T>* NRGA_getShapePoints_R4D();
	
	void NRGA_setCurrIterNo(long int currIterNo);
	void NRGA_simulate();
	
	void NRGA_Gravitational_Force2D_ARAP3(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D_NearPlane(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtreeT, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D_ARAP1(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D_ARAP2(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D_ARAP3(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Gravitational_Force3D_ARAP4(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
        void NRGA_Gravitational_Force3D_ARAP5(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, kdt_nanoFlann<T>* kdtree, kdt_nanoFlann<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_CurveSpaceGravitational_Force4D(point4D<T>* Tpts, point4D<T>* Rpts, point4D<T>* Tpts_prev, kdt_nanoFlann4D<T>* kdtree, kdt_nanoFlann4D<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);
	void NRGA_CurveSpaceGravitational2_Force4D(point4D<T>* Tpts, point4D<T>* Rpts, point4D<T>* Tpts_prev, kdt_nanoFlann4D<T>* kdtree, kdt_nanoFlann4D<T>* kdtreeR, int p_numPoints_T, int p_numPoints_R);

		
	void NRGA_Gravitational_Force2D(point3D<T>* Tpts, point3D<T>* Rpts, point3D<T>* Tpts_prev, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Spring_Force2D(point3D<T>* Tpts, point3D<T>* Tpts_prev, int p_numPoints_T, int p_numPoints_R);
	void NRGA_Spring_Force3D(point3D<T>* Tpts, point3D<T>* Tpts_prev, Eigen::Matrix<T, 3, 1>& centroid_template_curr, Eigen::Matrix<T, 3, 1>& centroid_template_prev,  int p_numPoints_T, int p_numPoints_R);
	void NRGA_sineDistort(point3D<T>* p, int numPts);
	void NRGA_calculateMeanDeviation(point3D<T>* Tpts, point3D<T>* Rpts, int numPtsTR);
	void NRGA_calculateRMSE(point3D<T>* Tpts, point3D<T>* Rpts, int numPtsTR, T meanDiv);
	void NRGA_save_as_flo(std::string save_path_flo, point3D<T>* Tpts_prev, point3D<T>* Tpts_curr, int N);
	void plotGPE_vs_Iteration(Gnuplot& gpe, std::vector<T>& gpe_values);


private:
	shape_point<T>* Tpt;
	point3D<T>* template_pts_prev3D;
	point3D<T>* template_points3D;
	point3D<T>* reference_points3D;
	point4D<T>* template_pts_prev4D;
	point4D<T>* template_points4D;
	point4D<T>* reference_points4D;
	long int currIterNo;
};

#endif
