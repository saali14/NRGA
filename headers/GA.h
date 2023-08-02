#ifndef _GA_
#define _GA_

/*
 * Gravitational Approach For PointSet Registration
 * @Author  : Sk Aziz Ali
 * @Date	: 5th May 2016
 * @Place	: DFKI, Kaiserslautern, Germany
*/

#pragma once
#include "parameters.h"
#include "gnuplot-iostream.h"
#include "KD_Tree.h"

template<typename T>
class GA
{
public:

	GA();
	GA(point3D<T>*, point3D<T>*, point3D<T>*);
	virtual ~GA();

	Gnuplot gpe, divg;
	std::vector<T> gpe_values, divg_values;


	void GA_setT(shape_point<T>* Tpt);
	shape_point<T>* GA_getT();
	void  GA_setParameters(parameters<T>* param);
	parameters<T>* GA_getParameters();

	void GA_initShapePoints_T(shape_point<T> &);
	void GA_initShapePoints_R(shape_point<T> &);
	void GA_initShapePoints_TDS(shape_point<T> &);
	void GA_initShapePoints_RDS(shape_point<T> &);

	//----------------------------------Rigid Body Cases(Latest Force Function)--------------------------------------------------------------------//
	void GA_gravitySimulation2D_latest(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_latest(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation2D_Efficient_latest(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_Efficient_latest(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation2D_NearFarPlane_latest(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_NearFarPlane_latest(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_subsampled_Efficient_latest(point3D<T>*, point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R, int p_numPoints_TDS, int p_numPoints_RDS);
	void GA_gravitySimulation3D_subsampled_Efficient_NearFarPlane_latest(point3D<T>*, point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R, int p_numPoints_TDS, int p_numPoints_RDS);
	//---------------------------------------------------------------------------------------------------------------------------------------------//



	//----------------------------------Rigid Body Cases--------------------------------------------------------------------------------------------//
	void GA_gravitySimulation2D(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation2D_Efficient(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_Efficient(point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation2D_NearFarPlane(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_NearFarPlane(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation2D_Efficient_NearFarPlane(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_Efficient_NearFarPlane(point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R);
	void GA_gravitySimulation3D_subsampled_Efficient(point3D<T>*, point3D<T>*, point3D<T>*, point3D<T>*, int p_numPoints_T, int p_numPoints_R, int p_numPoints_TDS, int p_numPoints_RDS);
	void GA_gravitySimulation3D_subsampled_Efficient_NearFarPlane(point3D<T>*, point3D<T>*, point3D<T>*, point3D<T>*, kdt_nanoFlann<T>*, int p_numPoints_T, int p_numPoints_R, int p_numPoints_TDS, int p_numPoints_RDS);
	//-----------------------------------------------------------------------------------------------------------------------------------------------//

	T    GA_divergenceX(T distX, T distY, T distZ, T massT, T massR);
	T    GA_divergenceY(T distX, T distY, T distZ, T massT, T massR);
	T    GA_divergenceZ(T distX, T distY, T distZ, T massT, T massR);
	T    GA_divergence_numerical_X(T distX, T distY, T distZ, T massT, T massR);
	T    GA_divergence_numerical_Y(T distX, T distY, T distZ, T massT, T massR);
	T    GA_divergence_numerical_Z(T distX, T distY, T distZ, T massT, T massR);
	
	
	T    GA_calculateMeanDeviation(point3D<T>*, point3D<T>*, long numPtsTR);
	T    GA_calculateRMSE(point3D<T>*, point3D<T>*, long numPtsTR, T meanDiv);
	void GA_normalizeFrame(point3D<T>*, point3D<T>*);
	void GA_denormalizeFrame(point3D<T>*, point3D<T>*);
	bool GA_isOutisdeBall(T, T);
	bool GA_isOutisdeBall3D(T, T, T);
	int  GA_getFrame(std::time_t startTime, std::time_t endTime);
	void GA_plotDivergence(Gnuplot&, std::vector<T>&);
	void GA_plotGPE_vs_Iteration(Gnuplot&, std::vector<T>&);
	void GA_saveTemplateReference(point3D<T>*, long, point3D<T>*, long);
	void GA_setCurrIterNo(int );

private:
	parameters<T>* param;
	shape_point<T>* Tpt;
	int currIterNo;
};

#endif // _GA_
