/*
* @author : Sk Aziz Ali
* @Date   : 21st July, 2016
* @Place  : Kolkata, India
*/

#include "NRGA.h"

template<class T> NRGA<T>::NRGA()
{
}

template<class T> void NRGA<T>::NRGA_setT(shape_point<T>* Tpt)
{
	this->Tpt = Tpt;
}

template<class T> shape_point<T>* NRGA<T>::NRGA_getT()
{
	return this->Tpt;
}

template<class T> void  NRGA<T>::NRGA_setParameters(parameters<T>* param)
{
	this->param = param;
}

template<class T> parameters<T>* NRGA<T>::NRGA_getParameters()
{
	return this->param;
}

template<class T> void  NRGA<T>::NRGA_setCurrIterNo(long int currIter)
{
	this->currIterNo = currIter;
}

template<class T> void NRGA<T>::NRGA_simulate()
{}


// initialize template points
template<class T> void NRGA<T>::NRGA_initShapePoints_T3D(shape_point<T> &template_pt)
{
    // initial template points
    this->template_pts_prev3D = (point3D<T>*)malloc(template_pt.pointSets3D.size() * sizeof(point3D<T>));
    this->template_points3D = (point3D<T>*)malloc(template_pt.pointSets3D.size() * sizeof(point3D<T>));

    for (int i = 0; i < template_pt.pointSets3D.size(); i++)
    {
	    this->template_points3D[i].pos.x = template_pt.pointSets3D.find(i)->second.pos.x;
	    this->template_points3D[i].pos.y = template_pt.pointSets3D.find(i)->second.pos.y;
	    this->template_points3D[i].pos.z = template_pt.pointSets3D.find(i)->second.pos.z;

	    this->template_points3D[i].v.vx = template_pt.pointSets3D.find(i)->second.v.vx;
	    this->template_points3D[i].v.vy = template_pt.pointSets3D.find(i)->second.v.vy;
	    this->template_points3D[i].v.vz = template_pt.pointSets3D.find(i)->second.v.vz;

	    this->template_points3D[i].a.ax = template_pt.pointSets3D.find(i)->second.a.ax;
	    this->template_points3D[i].a.ay = template_pt.pointSets3D.find(i)->second.a.ay;
	    this->template_points3D[i].a.az = template_pt.pointSets3D.find(i)->second.a.az;

	    this->template_points3D[i].color.R = template_pt.pointSets3D.find(i)->second.color.R;
	    this->template_points3D[i].color.G = template_pt.pointSets3D.find(i)->second.color.G;
	    this->template_points3D[i].color.B = template_pt.pointSets3D.find(i)->second.color.B;

	    this->template_points3D[i].m = template_pt.pointSets3D.find(i)->second.m;
    }

}

// initialize reference points
template<class T> void NRGA<T>::NRGA_initShapePoints_R3D(shape_point<T> &reference_pt)
{
    this->reference_points3D = (point3D<T>*)malloc(reference_pt.pointSets3D.size() * sizeof(point3D<T>));

    for (int i = 0; i < reference_pt.pointSets3D.size(); i++)
    {
	    this->reference_points3D[i].pos.x = reference_pt.pointSets3D.find(i)->second.pos.x;
	    this->reference_points3D[i].pos.y = reference_pt.pointSets3D.find(i)->second.pos.y;
	    this->reference_points3D[i].pos.z = reference_pt.pointSets3D.find(i)->second.pos.z;

	    this->reference_points3D[i].v.vx = reference_pt.pointSets3D.find(i)->second.v.vx;
	    this->reference_points3D[i].v.vy = reference_pt.pointSets3D.find(i)->second.v.vy;
	    this->reference_points3D[i].v.vz = reference_pt.pointSets3D.find(i)->second.v.vz;

	    this->reference_points3D[i].a.ax = reference_pt.pointSets3D.find(i)->second.a.ax;
	    this->reference_points3D[i].a.ay = reference_pt.pointSets3D.find(i)->second.a.ay;
	    this->reference_points3D[i].a.az = reference_pt.pointSets3D.find(i)->second.a.az;

	    this->reference_points3D[i].color.R = reference_pt.pointSets3D.find(i)->second.color.R;
	    this->reference_points3D[i].color.G = reference_pt.pointSets3D.find(i)->second.color.G;
	    this->reference_points3D[i].color.B = reference_pt.pointSets3D.find(i)->second.color.B;

	    this->reference_points3D[i].m = reference_pt.pointSets3D.find(i)->second.m;
    }

}

// initialize template points 
template<class T> void NRGA<T>::NRGA_initShapePoints_T4D(shape_point<T> &template_pt)
{
	// initial template points
	this->template_pts_prev4D = (point4D<T>*)malloc(template_pt.pointSets4D.size() * sizeof(point4D<T>));
	this->template_points4D = (point4D<T>*)malloc(template_pt.pointSets4D.size() * sizeof(point4D<T>));

	for (int i = 0; i < template_pt.pointSets4D.size(); i++)
	{
		this->template_points4D[i].pos.x = template_pt.pointSets4D.find(i)->second.pos.x;
		this->template_points4D[i].pos.y = template_pt.pointSets4D.find(i)->second.pos.y;
		this->template_points4D[i].pos.z = template_pt.pointSets4D.find(i)->second.pos.z;
		this->template_points4D[i].kappa = template_pt.pointSets4D.find(i)->second.kappa + (rand() / (RAND_MAX / 0.0001));

		this->template_points4D[i].v.vx = template_pt.pointSets4D.find(i)->second.v.vx;
		this->template_points4D[i].v.vy = template_pt.pointSets4D.find(i)->second.v.vy;
		this->template_points4D[i].v.vz = template_pt.pointSets4D.find(i)->second.v.vz;
		this->template_points4D[i].kappa_vel = template_pt.pointSets4D.find(i)->second.kappa_vel;

		this->template_points4D[i].a.ax = template_pt.pointSets4D.find(i)->second.a.ax;
		this->template_points4D[i].a.ay = template_pt.pointSets4D.find(i)->second.a.ay;
		this->template_points4D[i].a.az = template_pt.pointSets4D.find(i)->second.a.az;
		this->template_points4D[i].kappa_acc = template_pt.pointSets4D.find(i)->second.kappa_acc;

		this->template_points4D[i].color.R = template_pt.pointSets4D.find(i)->second.color.R;
		this->template_points4D[i].color.G = template_pt.pointSets4D.find(i)->second.color.G;
		this->template_points4D[i].color.B = template_pt.pointSets4D.find(i)->second.color.B;

		this->template_points4D[i].m = template_pt.pointSets4D.find(i)->second.m;
	}

}

// initialize reference points
template<class T> void NRGA<T>::NRGA_initShapePoints_R4D(shape_point<T> &reference_pt)
{
	this->reference_points4D = (point4D<T>*)malloc(reference_pt.pointSets4D.size() * sizeof(point4D<T>));

	for (int i = 0; i < reference_pt.pointSets4D.size(); i++)
	{
		this->reference_points4D[i].pos.x = reference_pt.pointSets4D.find(i)->second.pos.x;
		this->reference_points4D[i].pos.y = reference_pt.pointSets4D.find(i)->second.pos.y;
		this->reference_points4D[i].pos.z = reference_pt.pointSets4D.find(i)->second.pos.z;
		this->reference_points4D[i].kappa = reference_pt.pointSets4D.find(i)->second.kappa + (rand() / (RAND_MAX / 0.0001));
		
		this->reference_points4D[i].v.vx = reference_pt.pointSets4D.find(i)->second.v.vx;
		this->reference_points4D[i].v.vy = reference_pt.pointSets4D.find(i)->second.v.vy;
		this->reference_points4D[i].v.vz = reference_pt.pointSets4D.find(i)->second.v.vz;
		this->reference_points4D[i].kappa_vel = reference_pt.pointSets4D.find(i)->second.kappa_vel;
		
		this->reference_points4D[i].a.ax = reference_pt.pointSets4D.find(i)->second.a.ax;
		this->reference_points4D[i].a.ay = reference_pt.pointSets4D.find(i)->second.a.ay;
		this->reference_points4D[i].a.az = reference_pt.pointSets4D.find(i)->second.a.az;
		this->reference_points4D[i].kappa_acc = reference_pt.pointSets4D.find(i)->second.kappa_acc;
		
		this->reference_points4D[i].color.R = reference_pt.pointSets4D.find(i)->second.color.R;
		this->reference_points4D[i].color.G = reference_pt.pointSets4D.find(i)->second.color.G;
		this->reference_points4D[i].color.B = reference_pt.pointSets4D.find(i)->second.color.B;

		this->reference_points4D[i].m = reference_pt.pointSets4D.find(i)->second.m;
	}
}

template<class T> point4D<T>* NRGA<T>::NRGA_getShapePoints_T4DPrev()
{
  return this->template_pts_prev4D;
}

template<class T> point3D<T>* NRGA<T>::NRGA_getShapePoints_T3DPrev()
{
  return this->template_pts_prev3D;
}

template<class T> point3D<T>* NRGA<T>::NRGA_getShapePoints_T3D()
{
  return this->template_points3D;
}

template<class T> point3D<T>* NRGA<T>::NRGA_getShapePoints_R3D()
{
  return this->reference_points3D;
}

template<class T> point4D<T>* NRGA<T>::NRGA_getShapePoints_T4D()
{
  return this->template_points4D;
}

template<class T> point4D<T>* NRGA<T>::NRGA_getShapePoints_R4D()
{
  return this->reference_points4D;
}

// With Collective Motion
template<class T> void NRGA<T>::NRGA_Gravitational_Force3D_ARAP2(point3D<T>* template_points_curr, point3D<T>* reference_points, point3D<T>* template_pts_prev, kdt_nanoFlann<T>* KDTree_T, kdt_nanoFlann<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;
		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;
		template_pts_prev[t].m = template_points_curr[t].m;
	}

	this->NRGA_getT()->updateShapePoints3D(template_pts_prev, (long)p_numPoints_T);
	ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(3 * p_numPoints_T, 4);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 3);

		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;

			p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
				(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

			p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
				(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

			p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
				(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));
		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;
		p_accZ = -p_forceZ / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);
		
		template_points_curr[t].v.vx = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy = template_pts_prev[t].v.vy + p_velY;
		template_points_curr[t].v.vz = template_pts_prev[t].v.vz + p_velZ;

		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;

		T1(KDTree_T->k, 0) = template_points_curr[t].pos.x + (p_velX * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_points_curr[t].pos.y + (p_velY * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_points_curr[t].pos.z + (p_velZ * this->param->p_timeStep);

		/*
		Eigen::Matrix<T, 3, 1> avgVel(	(p_velX + template_points_curr[t].v.vx) / (KDTree_T->k + 1),
		(p_velY + template_points_curr[t].v.vy) / (KDTree_T->k + 1),
		(p_velZ + template_points_curr[t].v.vz) / (KDTree_T->k + 1)
		);
		
		T1(KDTree_T->k, 0) = template_points_curr[t].pos.x + ((avgVel(0, 0) /  avgVel.norm()) * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_points_curr[t].pos.y + ((avgVel(1, 0) /  avgVel.norm()) * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_points_curr[t].pos.z + ((avgVel(2, 0) /  avgVel.norm()) * this->param->p_timeStep);
		*/
		
		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			p_forceX = 0.0;
			p_forceY = 0.0;
			p_forceZ = 0.0;

			p_accX = 0.0;
			p_accY = 0.0;
			p_accZ = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;

				p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

				p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

				p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;

			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;

			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
		}

		this->NRGA_getT()->absOrientation(T0, T1, KDTree_T->k + 1, KDTree_T->k + 1, t, indexesNNT);

		/*
		Eigen::Matrix<T, 3, 1> cognitiveVel(((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x - template_points_curr[t].pos.x)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y - template_points_curr[t].pos.y)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z - template_points_curr[t].pos.z))
		);
		this->NRGA_getT()->updatePointVelocity3D(	template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vz + (cognitiveVel(2, 0) / cognitiveVel.norm()),
		t);
		this->NRGA_getT()->updatePointPosition(	this->param->p_timeStep * (template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vz + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		t);
		*/

		/*this->NRGA_getT()->updatePointVelocity3D( (avgVel(0, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		(avgVel(1, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		(avgVel(2, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		t);*/

		TransformationMap.block(3 * t, 0, 3, 4) = this->NRGA_getT()->getTransformation();
		outfile << this->NRGA_getT()->getTransformation();
		outfile << std::endl;

		//for (size_t _t = 0; _t < indexesNNT.size(); _t++)
		//{
		//this->NRGA_getT()->updatePointVelocity3D(	(avgVel(0, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		//										(avgVel(1, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		//										(avgVel(2, 0) / ((KDTree_T->k + 1) * avgVel.norm())),
		//										indexesNNT.at(_t));
		//}
	}

	// Produce Collective Motion by Modifying local velocitie
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_points_curr[t], indexesNNT, distNNT);

		Eigen::Matrix<T, 3, 1> avgVel( (template_points_curr[t].v.vx), (template_points_curr[t].v.vy), (template_points_curr[t].v.vz));

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
			avgVel(0,0) += template_points_curr[indexesNNT.at(i)].v.vx;
			avgVel(1,0) += template_points_curr[indexesNNT.at(i)].v.vy;
			avgVel(2,0) += template_points_curr[indexesNNT.at(i)].v.vz;
		}

		avgVel(0, 0) /= (KDTree_T->k + 1);
		avgVel(1, 0) /= (KDTree_T->k + 1);
		avgVel(2, 0) /= (KDTree_T->k + 1);

		this->NRGA_getT()->updatePointVelocity3D( (avgVel(0, 0) / avgVel.norm()), (avgVel(1,0)/avgVel.norm()), (avgVel(2,0)/avgVel.norm()), t);
	}

	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_points_curr[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
			this->NRGA_getT()->rotatePointXYZ(TransformationMap.block(3 * t, 0, 3, 3), indexesNNT.at(i));
			this->NRGA_getT()->translatePoint(TransformationMap.block(3 * t, 3, 3, 1), indexesNNT.at(i));
		}

		this->NRGA_getT()->rotatePointXYZ(TransformationMap.block(3 * t, 0, 3, 3), t);
		this->NRGA_getT()->translatePoint(TransformationMap.block(3 * t, 3, 3, 1), t);
	}


	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_points_curr[t].pos.x = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x;
		template_points_curr[t].pos.y = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y;
		template_points_curr[t].pos.z = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z;

		template_points_curr[t].v.vx = 0.95 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vx;
		template_points_curr[t].v.vy = 0.95 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vy;
		template_points_curr[t].v.vz = 0.95 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vz;
	}

	outfile.close();
	TransformationMap.resize(0, 0);

}

// With Collective Motion
template<class T> void NRGA<T>::NRGA_Gravitational_Force3D_ARAP3(point3D<T>* template_points_curr, point3D<T>* reference_points, point3D<T>* template_pts_prev, kdt_nanoFlann<T>* KDTree_T, kdt_nanoFlann<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
        T eps = this->param->p_epsilon;// * (exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
        T G   = this->param->p_G;//* (1.0 + exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;

		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;

		template_pts_prev[t].color.R = template_points_curr[t].color.R;
		template_pts_prev[t].color.G = template_points_curr[t].color.G;
		template_pts_prev[t].color.B = template_points_curr[t].color.B;

		template_pts_prev[t].m = template_points_curr[t].m;
	}

	//this->NRGA_getT()->updateShapePoints3D(template_pts_prev, (long)p_numPoints_T);
	//ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(3 * p_numPoints_T, 4);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;
                
                this->param->p_GPE = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);     
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 3);

		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;
			
			p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

			p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

			p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));
                                
                        this->param->p_GPE += -(this->param->p_G * template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / (std::sqrt((distX * distX) + (distY * distY) + (eps * eps)));

		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;
		p_accZ = -p_forceZ / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);

		template_points_curr[t].v.vx = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy = template_pts_prev[t].v.vy + p_velY;
		template_points_curr[t].v.vz = template_pts_prev[t].v.vz + p_velZ;

		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;

		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_points_curr[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_points_curr[t].v.vy * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_pts_prev[t].pos.z + (template_points_curr[t].v.vz * this->param->p_timeStep);


		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			p_forceX = 0.0;
			p_forceY = 0.0;
			p_forceZ = 0.0;

			p_accX = 0.0;
			p_accY = 0.0;
			p_accZ = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;

				p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

				p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

				p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));
                                        
                                this->param->p_GPE += -(this->param->p_G * template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / (std::sqrt((distX * distX) + (distY * distY) + (eps * eps)));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;

			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;

			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
		}

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, 3, 1> centroidT0;
		Eigen::Matrix<T, 3, 1> centroidT1;

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0) = (T)1.0f; S(1) = (T)1.0f;
		S(2) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - (Sk(0) * Rk * centroidT0));

		TransformationMap.block(3 * t, 0, 3, 3) = Rk;
		TransformationMap.block(3 * t, 3, 3, 1) = Tk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMapFS = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3 * p_numPoints_T, 5);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> frequency = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(p_numPoints_T, 1);
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
		    TransformationMapFS.block(3 * indexesNNT.at(i), 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		    TransformationMapFS(3 * indexesNNT.at(i) + 0, 4) +=  template_points_curr[t].v.vx;
		    TransformationMapFS(3 * indexesNNT.at(i) + 1, 4) +=  template_points_curr[t].v.vy;
		    TransformationMapFS(3 * indexesNNT.at(i) + 2, 4) +=  template_points_curr[t].v.vz;
		    frequency(indexesNNT.at(i)) += 1;
		}
		
		TransformationMapFS.block(3 * t, 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		frequency(t) += 1;
	}
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tpt_t(3,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y, template_pts_prev[t].pos.z;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 0, 3, 3)) * Tpt_t).transpose()				
					      + 
					     ((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 3, 3, 1)).transpose();
		
		template_points_curr[t].v.vx =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 0, 4));
		template_points_curr[t].v.vy =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 1, 4));
		template_points_curr[t].v.vz =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 2, 4));

		template_points_curr[t].pos.x = g_t(0);
		template_points_curr[t].pos.y = g_t(1);
		template_points_curr[t].pos.z = g_t(2);
		
		//template_points_curr[t].pos.x += template_points_curr[t].v.vx * this->param->p_timeStep;
		//template_points_curr[t].pos.y += template_points_curr[t].v.vy * this->param->p_timeStep;
		//template_points_curr[t].pos.z += template_points_curr[t].v.vz * this->param->p_timeStep;
		
		
// 		double f = std::sqrt( (template_points_curr[t].v.vx) * (template_points_curr[t].v.vx) 
// 				    + (template_points_curr[t].v.vy) * (template_points_curr[t].v.vy)  
// 				    + (template_points_curr[t].v.vz) * (template_points_curr[t].v.vz)
// 				      );
// 		double a = f/6.0f;//invert and group
// // 		double a =(1 - (f))/0.3;//invert and group
// 		int X = std::floor(a);	//this is the integer part
// 		int Y = std::floor(255*(a-X)); //fractional part from 0 to 255
// 		Eigen::VectorXi color_(3);
// 		
// 		if(a >= 0.0 && a < 0.25)
// 		  color_ << 0, 0, 128  +  std::floor(a * 255);
// 		else if(a >= 0.25 && a < 0.5)
// 		  color_ << 0, 128 + std::floor(a * 255), 0;
// 		else if(a >= 0.5 && a < 0.75)
// 		  color_ << std::floor(a * 255), 255, 0;
// 		else if(a >= 0.75 && a < 0.95)
// 		  color_ << std::floor(a * 255), 128, 0;
// 		else if( a >= 0.95)
// 		  color_ << 255, 0, 0;
// 		
// 		template_points_curr[t].color.R = color_(0);
// 		template_points_curr[t].color.G = color_(1);
// 		template_points_curr[t].color.B = color_(2);

		template_points_curr[t].v.vx *= 0.8;
		template_points_curr[t].v.vy *= 0.8;
		template_points_curr[t].v.vz *= 0.8;		
	}
	
	//outfile.close();
	TransformationMap.resize(0, 0);
	TransformationMapFS.resize(0, 0);
        
        this->GPE_values.push_back(this->param->p_GPE);
        this->plotGPE_vs_Iteration(this->gp, this->GPE_values);
}




// -----------------------------------------WITH COLLECTIVE MOTION--------------------------------------------------------
template<class T> void NRGA<T>::NRGA_Gravitational_Force2D_ARAP3(point3D<T>* template_points_curr, point3D<T>* reference_points, point3D<T>* template_pts_prev, kdt_nanoFlann<T>* KDTree_T, kdt_nanoFlann<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
        T eps = this->param->p_epsilon * (exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
        T G   = this->param->p_G;//* (1.0 + exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;

		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;

		template_pts_prev[t].m = template_points_curr[t].m;
	}

	//this->NRGA_getT()->updateShapePoints3D(template_pts_prev, (long)p_numPoints_T);
	//ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(2 * p_numPoints_T, 3);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);     
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 2);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 2);

		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			
			p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY)  + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (eps * eps))));

			p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY)  + (this->param->p_epsilon * this->param->p_epsilon))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (eps * eps))));
		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);;

		template_points_curr[t].v.vx = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy = template_pts_prev[t].v.vy + p_velY;

		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;

		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_points_curr[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_points_curr[t].v.vy * this->param->p_timeStep);


		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			p_forceX = 0.0;
			p_forceY = 0.0;

			p_accX = 0.0;
			p_accY = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;

				p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY)  + (eps * eps))));

				p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (eps * eps))));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;

			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;

			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
		}

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(2, 2);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(2, 2);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(2, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(2, 2), W(2, 2);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, 2, 1> centroidT0;
		Eigen::Matrix<T, 2, 1> centroidT1;

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0) = (T)1.0f;
		S(1) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - (Sk(0) * Rk * centroidT0));

		TransformationMap.block(2 * t, 0, 2, 2) = Rk;
		TransformationMap.block(2 * t, 2, 2, 1) = Tk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMapFS = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2 * p_numPoints_T, 4);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> frequency = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(p_numPoints_T, 1);
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
		    TransformationMapFS.block(2 * indexesNNT.at(i), 0, 2, 3) +=  TransformationMap.block(2*t, 0, 2, 3);
		    TransformationMapFS(2 * indexesNNT.at(i) + 0, 3) +=  template_points_curr[t].v.vx;
		    TransformationMapFS(2 * indexesNNT.at(i) + 1, 3) +=  template_points_curr[t].v.vy;
		    TransformationMapFS(2 * indexesNNT.at(i) + 2, 3) +=  template_points_curr[t].v.vz;
		    frequency(indexesNNT.at(i)) += 1;
		}
		
		TransformationMapFS.block(2 * t, 0, 2, 3) +=  TransformationMap.block(2*t, 0, 2, 3);
		frequency(t) += 1;
	}
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tpt_t(2,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (((1.0f/frequency(t)) * TransformationMapFS.block(2 * t, 0, 2, 2)) * Tpt_t).transpose()				
					      + 
					     ((1.0f/frequency(t)) * TransformationMapFS.block(2 * t, 2, 2, 1)).transpose();
		
		template_points_curr[t].v.vx =  ((1.0f/frequency(t)) * TransformationMapFS(2 * t + 0, 3));// + ((g_t(0, 0) - template_pts_prev[t].pos.x))/10.0f;
		template_points_curr[t].v.vy =  ((1.0f/frequency(t)) * TransformationMapFS(2 * t + 1, 3));// + ((g_t(0, 1) - template_pts_prev[t].pos.y))/10.0f;

		template_points_curr[t].pos.x = g_t(0, 0); //template_pts_prev[t].pos.x + template_points_curr[t].v.vx * this->param->p_timeStep;
		template_points_curr[t].pos.y = g_t(0, 1); //template_pts_prev[t].pos.y + template_points_curr[t].v.vy * this->param->p_timeStep;
		
		
		template_points_curr[t].v.vx *= 0.9;
		template_points_curr[t].v.vy *= 0.9;
	}

	//outfile.close();
	TransformationMap.resize(0, 0);
	TransformationMapFS.resize(0, 0);

}




// -------------------------------------WITH COHERENT COLLECTIVE MOTION ------------------------------------------------
template<class T> void NRGA<T>::NRGA_Gravitational_Force3D_ARAP4(point3D<T>* template_points_curr, point3D<T>* reference_points, point3D<T>* template_pts_prev, kdt_nanoFlann<T>* KDTree_T, kdt_nanoFlann<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;
		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;
		template_pts_prev[t].m = template_points_curr[t].m;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(3 * p_numPoints_T, 4);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;
		
		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 3);

		// Attractive Force by Gravitation 
		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;

			p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

			p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

			p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));	
		}
		
		// Repulsive Force module by Viscek Model 
		for(int s = 0; s < KDTree_T->k; s++)
		{
		  
			  distX =  template_pts_prev[t].pos.x - template_pts_prev[indexesNNT.at(s)].pos.x;
			  distY =  template_pts_prev[t].pos.y - template_pts_prev[indexesNNT.at(s)].pos.y;
			  distZ =  template_pts_prev[t].pos.z - template_pts_prev[indexesNNT.at(s)].pos.z;
			  
			  p_forceX +=  ((T)(-distX / (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) * 
					      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
			
			  p_forceY +=  ((T)(-distY / (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) *
					      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
			    
			  p_forceZ +=  ((T)(-distZ / (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) *
					      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
			  
			    /*
				  p_forceX +=  ((T)(-distX * 0.05)) * 
						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));
			
				  p_forceY +=  ((T)(-distY * 0.05)) *
						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));
				    
				  p_forceZ +=  ((T)(-distZ * 0.05)) *
						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));*/
		}
		
		p_accX = (-p_forceX )/ template_pts_prev[t].m;
		p_accY = (-p_forceY )/ template_pts_prev[t].m;
		p_accZ = (-p_forceZ )/ template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);

		template_pts_prev[t].v.vx = template_pts_prev[t].v.vx + p_velX;
		template_pts_prev[t].v.vy = template_pts_prev[t].v.vy + p_velY;
		template_pts_prev[t].v.vz = template_pts_prev[t].v.vz + p_velZ;

		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;

		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_pts_prev[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_pts_prev[t].v.vy * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_pts_prev[t].pos.z + (template_pts_prev[t].v.vz * this->param->p_timeStep);


		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			p_forceX = 0.0;
			p_forceY = 0.0;
			p_forceZ = 0.0;

			p_accX = 0.0;
			p_accY = 0.0;
			p_accZ = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;

				p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

				p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));

				p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))))) *
					(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (this->param->p_epsilon * this->param->p_epsilon))));
			}

			for(int s = 0; s < KDTree_T->k; s++)
			{
			
				 distX =  template_pts_prev[indexesNNT.at(_t)].pos.x - template_pts_prev[indexesNNT.at(s)].pos.x;
				 distY =  template_pts_prev[indexesNNT.at(_t)].pos.y - template_pts_prev[indexesNNT.at(s)].pos.y;
				 distZ =  template_pts_prev[indexesNNT.at(_t)].pos.z - template_pts_prev[indexesNNT.at(s)].pos.z;
				  
				 p_forceX +=  ((T)(-distX / (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) * 
					      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
			
				  p_forceY +=  ((T)(-distY/ (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) *
						      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
				    
				  p_forceZ +=  ((T)(-distZ/ (sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (this->param->p_epsilon * this->param->p_epsilon)))) *
						      (1.0 / (1.0 + exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/1.5) - 2))));
				  
// 				  p_forceX +=  ((T)(-distX * 0.05)) * 
// 						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));
// 			
// 				  p_forceY +=  ((T)(-distY * 0.05)) *
// 						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));
// 				    
// 				  p_forceZ +=  ((T)(-distZ * 0.05)) *
// 						      (1+  (exp( (sqrt((distX * distX) + (distY * distY) + (distZ * distZ))/0.2))));
				  
			}
			
			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;

			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;

			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
		}

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, 3, 1> centroidT0;
		Eigen::Matrix<T, 3, 1> centroidT1;

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0) = (T)1.0f; S(1) = (T)1.0f;
		S(2) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - (Sk(0) * Rk * centroidT0));

		TransformationMap.block(3 * t, 0, 3, 3) = Rk;
		TransformationMap.block(3 * t, 3, 3, 1) = Tk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);

		//std::cout << "reached here --->" << std::endl;
		//TransformationMap.block(3 * t, 0, 3, 4) = this->NRGA_getT()->getTransformation();
		//outfile << this->NRGA_getT()->getTransformation();
		//outfile << std::endl;
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMapFS = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3 * p_numPoints_T, 5);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> frequency = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(p_numPoints_T, 1);
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
		    TransformationMapFS.block(3 * indexesNNT.at(i), 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		    
		    TransformationMapFS(3 * indexesNNT.at(i) + 0, 4) +=  template_pts_prev[t].v.vx;
		    TransformationMapFS(3 * indexesNNT.at(i) + 1, 4) +=  template_pts_prev[t].v.vy;
		    TransformationMapFS(3 * indexesNNT.at(i) + 2, 4) +=  template_pts_prev[t].v.vz;
		    frequency(indexesNNT.at(i)) += 1;
		}
		
		TransformationMapFS(3 * t + 0, 4) +=  template_pts_prev[t].v.vx;
		TransformationMapFS(3 * t + 1, 4) +=  template_pts_prev[t].v.vy;
		TransformationMapFS(3 * t + 2, 4) +=  template_pts_prev[t].v.vz;
		TransformationMapFS.block(3 * t, 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		frequency(t) += 1;
	}
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tpt_t(3,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y, template_pts_prev[t].pos.z;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 0, 3, 3)) * Tpt_t).transpose()				
					      + 
					     ((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 3, 3, 1)).transpose();
		
		T v_0 = sqrt((template_pts_prev[t].v.vx * template_pts_prev[t].v.vx) + 
			     (template_pts_prev[t].v.vy * template_pts_prev[t].v.vy) + 
			     (template_pts_prev[t].v.vz * template_pts_prev[t].v.vz) 
			    );
		
		template_pts_prev[t].v.vx =  v_0 * ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 0, 4)/TransformationMapFS.block(3 * t, 4, 3, 1).norm());
		template_pts_prev[t].v.vy =  v_0 * ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 1, 4)/TransformationMapFS.block(3 * t, 4, 3, 1).norm());
		template_pts_prev[t].v.vz =  v_0 * ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 2, 4)/TransformationMapFS.block(3 * t, 4, 3, 1).norm());

// 		template_pts_prev[t].v.vx =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 0, 4));
//  		template_pts_prev[t].v.vy =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 1, 4));
//  		template_pts_prev[t].v.vz =  ((1.0f/frequency(t)) * TransformationMapFS(3 * t + 2, 4));

		template_points_curr[t].v.vx =  template_pts_prev[t].v.vx ;
		template_points_curr[t].v.vy =  template_pts_prev[t].v.vy ;
		template_points_curr[t].v.vz =  template_pts_prev[t].v.vz ;

//              template_points_curr[t].pos.x +=  (template_points_curr[t].v.vx * this->param->p_timeStep);
// 		template_points_curr[t].pos.y +=  (template_points_curr[t].v.vy * this->param->p_timeStep);
// 		template_points_curr[t].pos.z +=  (template_points_curr[t].v.vz * this->param->p_timeStep);
		
		template_points_curr[t].pos.x =  g_t(0, 0) ;//+ (template_points_curr[t].v.vx * this->param->p_timeStep);
		template_points_curr[t].pos.y =  g_t(0, 1) ;//+ (template_points_curr[t].v.vy * this->param->p_timeStep);
		template_points_curr[t].pos.z =  g_t(0, 2) ;//+ (template_points_curr[t].v.vz * this->param->p_timeStep);

		template_points_curr[t].v.vx *= 0.8;
		template_points_curr[t].v.vy *= 0.8;
		template_points_curr[t].v.vz *= 0.8;
	}
	
	//outfile.close();
	TransformationMap.resize(0, 0);
	TransformationMapFS.resize(0, 0);
}



// ------------------------------------WITHOUT COHERENT COLLECTIVE MOTION ----------------------------------------------
template<class T> void NRGA<T>::NRGA_Gravitational_Force3D_ARAP5(point3D<T>* template_points_curr, point3D<T>* reference_points, point3D<T>* template_pts_prev, kdt_nanoFlann<T>* KDTree_T, kdt_nanoFlann<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
T eps = this->param->p_epsilon;// * (exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
        T G   = this->param->p_G;//* (1.0 + exp(-(this->currIterNo/(this->param->p_max_iterations + 2))));
	
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;

		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;

		template_pts_prev[t].color.R = template_points_curr[t].color.R;
		template_pts_prev[t].color.G = template_points_curr[t].color.G;
		template_pts_prev[t].color.B = template_points_curr[t].color.B;

		template_pts_prev[t].m = template_points_curr[t].m;
	}

	//this->NRGA_getT()->updateShapePoints3D(template_pts_prev, (long)p_numPoints_T);
	//ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(3 * p_numPoints_T, 4);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);     
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 3);

		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;
			
			p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

			p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

			p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
				(G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / ((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;
		p_accZ = -p_forceZ / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);

		template_points_curr[t].v.vx = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy = template_pts_prev[t].v.vy + p_velY;
		template_points_curr[t].v.vz = template_pts_prev[t].v.vz + p_velZ;

		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;

		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_points_curr[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_points_curr[t].v.vy * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_pts_prev[t].pos.z + (template_points_curr[t].v.vz * this->param->p_timeStep);


		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			p_forceX = 0.0;
			p_forceY = 0.0;
			p_forceZ = 0.0;

			p_accX = 0.0;
			p_accY = 0.0;
			p_accZ = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;

				p_forceX += (T)(distX / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

				p_forceY += (T)(distY / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));

				p_forceZ += (T)(distZ / (sqrt(((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))))) *
					(G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) /
					((distX * distX) + (distY * distY) + (distZ * distZ) + (eps * eps))));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;

			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;

			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
		}

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, 3, 1> centroidT0;
		Eigen::Matrix<T, 3, 1> centroidT1;

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0) = (T)1.0f; S(1) = (T)1.0f;
		S(2) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - (Sk(0) * Rk * centroidT0));

		TransformationMap.block(3 * t, 0, 3, 3) = Rk;
		TransformationMap.block(3 * t, 3, 3, 1) = Tk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);
	}
	
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tpt_t(3,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y, template_pts_prev[t].pos.z;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (TransformationMap.block(3 * t, 0, 3, 3) * Tpt_t).transpose()	
					      + 
					     (TransformationMap.block(3 * t, 3, 3, 1)).transpose();
		
		template_points_curr[t].pos.x = g_t(0);
		template_points_curr[t].pos.y = g_t(1);
		template_points_curr[t].pos.z = g_t(2);
		
		template_points_curr[t].v.vx *= 0.8;
		template_points_curr[t].v.vy *= 0.8;
		template_points_curr[t].v.vz *= 0.8;		
	}
	
	//outfile.close();
	TransformationMap.resize(0, 0);

}

template<class T> void NRGA<T>::NRGA_CurveSpaceGravitational_Force4D(point4D<T>* template_points_curr, point4D<T>* reference_points, point4D<T>* template_pts_prev, kdt_nanoFlann4D<T>* KDTree_T, kdt_nanoFlann4D<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
	///@COMMENTS:: First define the new General Relativistic Field Equation(Einstein's Field Equation)
	///            and reformulation of post Newtonian Gravitational Force function using 
	///		--> G as a function of "Curvature Tensor/ 8 * Pi * (Constant)" 
	///		--> G as a function of "Torsion Tensor "
	///		--> Template masses a "Radial Basis Function of Absolute Curvature Difference"
	///
	///@NEW_Gratitational_Force_Formulation-->
	///@UNDER_VISCOUS_MEDIUM--> The space is assumed to be a viscous medium 
	///@PATCH_WISE_KABSCH--> Allign the X(patch = p, time =t) to X(patch = p, time =t+1)
	///@SAME_POINT_IN_MULTIPLE_PATCHES--> Same template point may belong to multiple patches
	///@REGULARIZED_GOAL_POSITION--> We use a Fast Summation from LAttice based shape matciing to 
	///				 uniformly apply one Transformation per template point
	///@COLLECTIVE_MOTION--> Updates the template points' velocity avoiding crossover and 
	///			 maintaing coherent motion (inspired by Self-Propelled Particle System)
  
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;
		template_pts_prev[t].kappa = template_points_curr[t].kappa;
		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;
		template_pts_prev[t].kappa_vel = template_points_curr[t].kappa_vel;
		template_pts_prev[t].m = template_points_curr[t].m;
	}

	this->NRGA_getT()->updateShapePoints4D(template_pts_prev, (long)p_numPoints_T);
	//ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(4 * p_numPoints_T, 5);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;
	T distK = 0.0;
	T r_ij	= 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;
		T p_forceK = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;
		T p_accK = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;
		T p_velK = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 4);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 4);


		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;
			distK = template_pts_prev[t].kappa - reference_points[indexesNNR.at(r)].kappa;
			
			if(reference_points[indexesNNR.at(r)].kappa < 0)
			{
			  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  - (distK * distK));
			}
			else if(reference_points[indexesNNR.at(r)].kappa >= 0)
			{
			  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (distK * distK));
			}
			else
			{
			  ; //r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ));
			}

			T k_facNumerator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/2.0f);
			T k_facDenomenator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/4.0f);
			
			p_forceX += (T)((template_pts_prev[t].pos.x * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.x) / 
					  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					  (((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceY += (T)((template_pts_prev[t].pos.y * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.y) / 
					  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceZ += (T)((template_pts_prev[t].pos.z * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.z) /
					sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceK += (T)((template_pts_prev[t].kappa * k_facNumerator) - reference_points[indexesNNR.at(r)].kappa) /
					sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));
		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;
		p_accZ = -p_forceZ / template_pts_prev[t].m;
		p_accK = -p_forceK / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);
		p_velK = (p_accK * this->param->p_timeStep);
		
		template_points_curr[t].v.vx 	  = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy 	  = template_pts_prev[t].v.vy + p_velY;
		template_points_curr[t].v.vz 	  = template_pts_prev[t].v.vz + p_velZ;
		template_points_curr[t].kappa_vel = template_pts_prev[t].kappa_vel + p_velK;
		
		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;
		T0(KDTree_T->k, 3) = template_pts_prev[t].kappa;
		
		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_points_curr[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_points_curr[t].v.vy * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_pts_prev[t].pos.z + (template_points_curr[t].v.vz * this->param->p_timeStep);
		T1(KDTree_T->k, 3) = template_pts_prev[t].kappa + (template_points_curr[t].kappa_vel * this->param->p_timeStep);

		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			//calculate force on each point using Gravitational Force field 
			T p_forceX = 0.0;
			T p_forceY = 0.0;
			T p_forceZ = 0.0;
			T p_forceK = 0.0;

			T p_accX = 0.0;
			T p_accY = 0.0;
			T p_accZ = 0.0;
			T p_accK = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;
				distK = template_pts_prev[indexesNNT.at(_t)].kappa - reference_points[indexesNNR.at(r)].kappa;
				
				if(reference_points[indexesNNR.at(r)].kappa < 0)
				{
				  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  - (distK * distK));
				}
				else if(reference_points[indexesNNR.at(r)].kappa >= 0)
				{
				  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (distK * distK));
				}
				else
				{
				  ;//r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ));
				}
				
				T k_facNumerator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/2.0f);
				T k_facDenomenator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/4.0f);
				
				p_forceX += (T)((template_pts_prev[indexesNNT.at(_t)].pos.x * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.x) / 
						  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						  (((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceY += (T)((template_pts_prev[indexesNNT.at(_t)].pos.y * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.y) / 
						  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceZ += (T)((template_pts_prev[indexesNNT.at(_t)].pos.z * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.z) /
						sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceK += (T)((template_pts_prev[indexesNNT.at(_t)].kappa * k_facNumerator) - reference_points[indexesNNR.at(r)].kappa) /
						sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;
			p_accK = -p_forceK / template_pts_prev[indexesNNT.at(_t)].m;
			
			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;
			T0(_t, 3) = template_pts_prev[indexesNNT.at(_t)].kappa;
			
			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 3) = template_pts_prev[indexesNNT.at(_t)].kappa + (p_accK * this->param->p_timeStep * this->param->p_timeStep);
		}
		
		//this->NRGA_getT()->absOrientation(T0, T1, KDTree_T->k + 1, KDTree_T->k + 1, t, indexesNNT);

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk(T0.cols(), 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk(T0.cols(), 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(T0.cols(), T0.cols()), W(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT0(T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT1(T1.cols());

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0,0) = (T)1.0f; S(1,0) = (T)1.0f;S(2,0) = (T)1.0f;
		S(3,0) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - Rk * centroidT0);

		TransformationMap.block(4*t, 0, 4, 4) = Rk;
		TransformationMap.block(4*t, 4, 4, 1) = Tk;
		//this->scale = Sk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);
	
		/*
		Eigen::Matrix<T, 3, 1> cognitiveVel(((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x - template_points_curr[t].pos.x)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y - template_points_curr[t].pos.y)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z - template_points_curr[t].pos.z))
		);
		this->NRGA_getT()->updatePointVelocity3D(	template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vz + (cognitiveVel(2, 0) / cognitiveVel.norm()),
		t);
		this->NRGA_getT()->updatePointPosition(	this->param->p_timeStep * (template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vz + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		t);
		*/
		//std::cout << "Reached Here -->" << std::endl;
		//TransformationMap.block(4 * t, 0, 4, 5) = this->NRGA_getT()->getTransformation();
		//outfile << this->NRGA_getT()->getTransformation();
		//outfile << std::endl;
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMapFS = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(4 * p_numPoints_T, 5);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> frequency = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(p_numPoints_T, 1);
	
	
	// IMPORTANT NOTE:: One doesnot know from howmany NETS or LATTICE_BLOCKS it belongs to.
	// CALCULATION of FREQUENCY HISTOGRAM  of every template point
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
		    TransformationMapFS.block(4 * indexesNNT.at(i), 0, 4, 5) +=  TransformationMap.block(4*t, 0, 4, 5);
		    frequency(indexesNNT.at(i)) += 1;
		}
		
		TransformationMapFS.block(4 * t, 0, 4, 5) +=  TransformationMap.block(4*t, 0, 4, 5);
		frequency(t) += 1;
	}
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Tpt_t(4,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y, template_pts_prev[t].pos.z, template_pts_prev[t].kappa;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (((1.0f/frequency(t)) * TransformationMapFS.block(4 * t, 0, 4, 4)) * Tpt_t).transpose()				
					      + 
					     ((1.0f/frequency(t)) * TransformationMapFS.block(4 * t, 4, 4, 1)).transpose();
		
// 		template_points_curr[t].v.vx +=  ((g_t(0, 0) - template_pts_prev[t].pos.x))/10.0f;
// 		template_points_curr[t].v.vy +=  ((g_t(0, 1) - template_pts_prev[t].pos.y))/10.0f;
// 		template_points_curr[t].v.vz +=  ((g_t(0, 2) - template_pts_prev[t].pos.z))/10.0f;
// 		template_points_curr[t].kappa_vel +=  ((g_t(0, 3) - template_pts_prev[t].kappa))/5.0f;

		template_points_curr[t].v.vx = (1.1 * template_points_curr[t].v.vx) + (g_t(0, 0) - template_pts_prev[t].pos.x)/20.0f;
		template_points_curr[t].v.vy = (1.1 * template_points_curr[t].v.vy) + (g_t(0, 1) - template_pts_prev[t].pos.y)/20.0f;
		template_points_curr[t].v.vz = (1.1 * template_points_curr[t].v.vz) + (g_t(0, 2) - template_pts_prev[t].pos.z)/20.0f;
		template_points_curr[t].kappa_vel = (1.1 * template_points_curr[t].kappa_vel) +  (g_t(0, 3) - template_pts_prev[t].kappa) /20.0f;
		
		template_points_curr[t].pos.x = template_pts_prev[t].pos.x + template_points_curr[t].v.vx * this->param->p_timeStep;
		template_points_curr[t].pos.y = template_pts_prev[t].pos.y + template_points_curr[t].v.vy * this->param->p_timeStep;
		template_points_curr[t].pos.z = template_pts_prev[t].pos.z + template_points_curr[t].v.vz * this->param->p_timeStep;
		template_points_curr[t].kappa = template_pts_prev[t].kappa + template_points_curr[t].kappa_vel * this->param->p_timeStep;
		
		template_points_curr[t].v.vx *= 0.1;
		template_points_curr[t].v.vy *= 0.1;
		template_points_curr[t].v.vz *= 0.1;
		//template_points_curr[t].kappa_vel *= 0.1;

// 		this->NRGA_getT()->updatePointVelocity3D(template_points_curr[t].v.vx, template_points_curr[t].v.vy, template_points_curr[t].v.vz, t);
// 		this->NRGA_getT()->updatePointPosition(template_points_curr[t].pos.x, template_points_curr[t].pos.y,template_points_curr[t].pos.z,t);
	}

// 	for (int t = 0; t < p_numPoints_T; t++)
// 	{
// 		template_points_curr[t].pos.x = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x;
// 		template_points_curr[t].pos.y = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y;
// 		template_points_curr[t].pos.z = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z;
// 
// 		template_points_curr[t].v.vx = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vx;
// 		template_points_curr[t].v.vy = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vy;
// 		template_points_curr[t].v.vz = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vz;
// 	}

	//outfile.close();
	TransformationMap.resize(0, 0);
	TransformationMapFS.resize(0, 0);
}

template<class T> void NRGA<T>::NRGA_CurveSpaceGravitational2_Force4D(point4D<T>* template_points_curr, point4D<T>* reference_points, point4D<T>* template_pts_prev, kdt_nanoFlann4D<T>* KDTree_T, kdt_nanoFlann4D<T>* KDTree_R, int p_numPoints_T, int p_numPoints_R)
{
	///@COMMENTS:: First define the new General Relativistic Field Equation(Einstein's Field Equation)
	///            and reformulation of post Newtonian Gravitational Force function using 
	///		--> G as a function of "Curvature Tensor/ 8 * Pi * (Constant)" 
	///		--> G as a function of "Torsion Tensor "
	///		--> Template masses a "Radial Basis Function of Absolute Curvature Difference"
	///
	///@NEW_Gratitational_Force_Formulation-->
	///@UNDER_VISCOUS_MEDIUM--> The space is assumed to be a viscous medium 
	///@PATCH_WISE_KABSCH--> Allign the X(patch = p, time =t) to X(patch = p, time =t+1)
	///@SAME_POINT_IN_MULTIPLE_PATCHES--> Same template point may belong to multiple patches
	///@REGULARIZED_GOAL_POSITION--> We use a Fast Summation from LAttice based shape matciing to 
	///				 uniformly apply one Transformation per template point
	///@COLLECTIVE_MOTION--> Updates the template points' velocity avoiding crossover and 
	///			 maintaing coherent motion (inspired by Self-Propelled Particle System)
  
	for (int t = 0; t < p_numPoints_T; t++)
	{
		template_pts_prev[t].pos.x = template_points_curr[t].pos.x;
		template_pts_prev[t].pos.y = template_points_curr[t].pos.y;
		template_pts_prev[t].pos.z = template_points_curr[t].pos.z;
		template_pts_prev[t].kappa = template_points_curr[t].kappa;
		template_pts_prev[t].v.vx = template_points_curr[t].v.vx;
		template_pts_prev[t].v.vy = template_points_curr[t].v.vy;
		template_pts_prev[t].v.vz = template_points_curr[t].v.vz;
		template_pts_prev[t].kappa_vel = template_points_curr[t].kappa_vel;
		template_pts_prev[t].m = template_points_curr[t].m;
	}

	this->NRGA_getT()->updateShapePoints4D(template_pts_prev, (long)p_numPoints_T);
	//ofstream outfile("iterTransformation_" + std::to_string(this->currIterNo) + ".txt", ios::app);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMap(3 * p_numPoints_T, 4);

	T distX = 0.0;
	T distY = 0.0;
	T distZ = 0.0;
	T distK = 0.0;
	T r_ij	= 0.0;

	for (int t = 0; t < p_numPoints_T; t++)
	{
		//calculate force on each point using Gravitational Force field 
		T p_forceX = 0.0;
		T p_forceY = 0.0;
		T p_forceZ = 0.0;
		T p_forceK = 0.0;

		T p_accX = 0.0;
		T p_accY = 0.0;
		T p_accZ = 0.0;
		T p_accK = 0.0;

		T p_velX = 0.0;
		T p_velY = 0.0;
		T p_velZ = 0.0;
		T p_velK = 0.0;

		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		std::vector<size_t> indexesNNR(KDTree_R->k);
		std::vector<T> distNNR(KDTree_R->k);
		KDTree_R->getKNNS(template_pts_prev[t], indexesNNR, distNNR);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T0(KDTree_T->k + 1, 3);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T1(KDTree_T->k + 1, 3);


		for (int r = 0; r < KDTree_R->k; ++r)
		{
			distX = template_pts_prev[t].pos.x - reference_points[indexesNNR.at(r)].pos.x;
			distY = template_pts_prev[t].pos.y - reference_points[indexesNNR.at(r)].pos.y;
			distZ = template_pts_prev[t].pos.z - reference_points[indexesNNR.at(r)].pos.z;
			distK = template_pts_prev[t].kappa - reference_points[indexesNNR.at(r)].kappa;
			
			if(reference_points[indexesNNR.at(r)].kappa < 0)
			{
			  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  - (distK * distK));
			}
			else if(reference_points[indexesNNR.at(r)].kappa >= 0)
			{
			  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (distK * distK));
			}
			else
			{
			  ; //r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ));
			}

			T k_facNumerator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/2.0f);
			T k_facDenomenator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/4.0f);
			
			p_forceX += (T)((template_pts_prev[t].pos.x * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.x) / 
					  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					  (((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceY += (T)((template_pts_prev[t].pos.y * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.y) / 
					  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceZ += (T)((template_pts_prev[t].pos.z * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.z) /
					sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

			p_forceK += (T)((template_pts_prev[t].kappa * k_facNumerator) - reference_points[indexesNNR.at(r)].kappa) /
					sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
					(this->param->p_G * ((template_pts_prev[t].m * reference_points[indexesNNR.at(r)].m) / 
					(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));
		}

		p_accX = -p_forceX / template_pts_prev[t].m;
		p_accY = -p_forceY / template_pts_prev[t].m;
		p_accZ = -p_forceZ / template_pts_prev[t].m;
		p_accK = -p_forceK / template_pts_prev[t].m;

		p_velX = (p_accX * this->param->p_timeStep);
		p_velY = (p_accY * this->param->p_timeStep);
		p_velZ = (p_accZ * this->param->p_timeStep);
		p_velK = (p_accK * this->param->p_timeStep);
		
		template_points_curr[t].v.vx 	  = template_pts_prev[t].v.vx + p_velX;
		template_points_curr[t].v.vy 	  = template_pts_prev[t].v.vy + p_velY;
		template_points_curr[t].v.vz 	  = template_pts_prev[t].v.vz + p_velZ;
		template_points_curr[t].kappa_vel = template_pts_prev[t].kappa_vel + p_velK;
		
		T0(KDTree_T->k, 0) = template_pts_prev[t].pos.x;
		T0(KDTree_T->k, 1) = template_pts_prev[t].pos.y;
		T0(KDTree_T->k, 2) = template_pts_prev[t].pos.z;
		//T0(KDTree_T->k, 3) = template_pts_prev[t].kappa;
		
		T1(KDTree_T->k, 0) = template_pts_prev[t].pos.x + (template_points_curr[t].v.vx * this->param->p_timeStep);
		T1(KDTree_T->k, 1) = template_pts_prev[t].pos.y + (template_points_curr[t].v.vy * this->param->p_timeStep);
		T1(KDTree_T->k, 2) = template_pts_prev[t].pos.z + (template_points_curr[t].v.vz * this->param->p_timeStep);
		//T1(KDTree_T->k, 3) = template_pts_prev[t].kappa + (template_points_curr[t].kappa_vel * this->param->p_timeStep);

		for (int _t = 0; _t < KDTree_T->k; _t++)
		{
			//calculate force on each point using Gravitational Force field 
			T p_forceX = 0.0;
			T p_forceY = 0.0;
			T p_forceZ = 0.0;
			T p_forceK = 0.0;

			T p_accX = 0.0;
			T p_accY = 0.0;
			T p_accZ = 0.0;
			T p_accK = 0.0;

			for (int r = 0; r < KDTree_R->k; ++r)
			{
				distX = template_pts_prev[indexesNNT.at(_t)].pos.x - reference_points[indexesNNR.at(r)].pos.x;
				distY = template_pts_prev[indexesNNT.at(_t)].pos.y - reference_points[indexesNNR.at(r)].pos.y;
				distZ = template_pts_prev[indexesNNT.at(_t)].pos.z - reference_points[indexesNNR.at(r)].pos.z;
				distK = template_pts_prev[indexesNNT.at(_t)].kappa - reference_points[indexesNNR.at(r)].kappa;
				
				if(reference_points[indexesNNR.at(r)].kappa < 0)
				{
				  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  - (distK * distK));
				}
				else if(reference_points[indexesNNR.at(r)].kappa >= 0)
				{
				  r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ)  + (distK * distK));
				}
				else
				{
				  ;//r_ij = sqrt((distX * distX) + (distY * distY) + (distZ * distZ));
				}
				
				T k_facNumerator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/2.0f);
				T k_facDenomenator = (T)(1.0f - (reference_points[indexesNNR.at(r)].kappa * r_ij * r_ij)/4.0f);
				
				p_forceX += (T)((template_pts_prev[indexesNNT.at(_t)].pos.x * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.x) / 
						  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						  (((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceY += (T)((template_pts_prev[indexesNNT.at(_t)].pos.y * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.y) / 
						  sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceZ += (T)((template_pts_prev[indexesNNT.at(_t)].pos.z * k_facNumerator) - reference_points[indexesNNR.at(r)].pos.z) /
						sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));

				p_forceK += (T)((template_pts_prev[indexesNNT.at(_t)].kappa * k_facNumerator) - reference_points[indexesNNR.at(r)].kappa) /
						sqrt((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) *
						(this->param->p_G * ((template_pts_prev[indexesNNT.at(_t)].m * reference_points[indexesNNR.at(r)].m) / 
						(((r_ij * r_ij) + (this->param->p_epsilon * this->param->p_epsilon)) * (k_facDenomenator * k_facDenomenator))));
			}

			p_accX = -p_forceX / template_pts_prev[indexesNNT.at(_t)].m;
			p_accY = -p_forceY / template_pts_prev[indexesNNT.at(_t)].m;
			p_accZ = -p_forceZ / template_pts_prev[indexesNNT.at(_t)].m;
			p_accK = -p_forceK / template_pts_prev[indexesNNT.at(_t)].m;
			
			T0(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x;
			T0(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y;
			T0(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z;
			//T0(_t, 3) = template_pts_prev[indexesNNT.at(_t)].kappa;
			
			T1(_t, 0) = template_pts_prev[indexesNNT.at(_t)].pos.x + (p_accX * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 1) = template_pts_prev[indexesNNT.at(_t)].pos.y + (p_accY * this->param->p_timeStep * this->param->p_timeStep);
			T1(_t, 2) = template_pts_prev[indexesNNT.at(_t)].pos.z + (p_accZ * this->param->p_timeStep * this->param->p_timeStep);
			//T1(_t, 3) = template_pts_prev[indexesNNT.at(_t)].kappa + (p_accK * this->param->p_timeStep * this->param->p_timeStep);
		}
		
		//this->NRGA_getT()->absOrientation(T0, T1, KDTree_T->k + 1, KDTree_T->k + 1, t, indexesNNT);

		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk(T0.cols(), 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk(T0.cols(), 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(T0.cols(), T0.cols()), W(T0.cols(), T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
		Eigen::Matrix<T, Eigen::Dynamic, 1> S;
		Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT0(T0.cols());
		Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT1(T1.cols());

		// Compute Centroid
		centroidT0 = T0.colwise().mean();
		centroidT1 = T1.colwise().mean();

		// Shift points by centroids 
		T0.rowwise() -= centroidT0.transpose();
		T1.rowwise() -= centroidT1.transpose();
		
		// Compute linear deformation Matrix -> Ak
		Ak = (T0.transpose() * T1) * (T0.transpose() * T0).inverse();

		// Compute Orientation -> Rk
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(T0.transpose() * T1, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();
		
		WVT = W * VT;
		S(0,0) = (T)1.0f; S(1,0) = (T)1.0f;
		S(2,0) = (T)(WVT.determinant() <= 0) ? -1.0f : 1.0f; 
		Rk = W * S.asDiagonal() * VT;


		// Compute Scale -> Sk
		T sT0 = (T0.array() * T0.array()).sum();
		T sT1 = (T1.array() * T1.array()).sum();
		Sk(0) = (T)1.0f; // std::sqrt(sT1 / sT0); Sk(1) = std::sqrt(sT1 / sT0); Sk(2) = std::sqrt(sT1 / sT0);

		// Compute Translation -> Tk
		Tk = (centroidT1 - Rk * centroidT0);

		TransformationMap.block(3*t, 0, 3, 3) = Rk;
		TransformationMap.block(3*t, 3, 3, 1) = Tk;
		//this->scale = Sk;
		
		Ak.resize(0, 0);
		Sk.resize(0, 0);
		Rk.resize(0, 0);
		Tk.resize(0, 0);
		VT.resize(0, 0);
		W.resize(0, 0);
		WVT.resize(0, 0);
	
		/*
		Eigen::Matrix<T, 3, 1> cognitiveVel(((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x - template_points_curr[t].pos.x)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y - template_points_curr[t].pos.y)),
		((this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z - template_points_curr[t].pos.z))
		);
		this->NRGA_getT()->updatePointVelocity3D(	template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm()),
		template_points_curr[t].v.vz + (cognitiveVel(2, 0) / cognitiveVel.norm()),
		t);
		this->NRGA_getT()->updatePointPosition(	this->param->p_timeStep * (template_points_curr[t].v.vx + (cognitiveVel(0, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vy + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		this->param->p_timeStep * (template_points_curr[t].v.vz + (cognitiveVel(1, 0) / cognitiveVel.norm())),
		t);
		*/
		//std::cout << "Reached Here -->" << std::endl;
		//TransformationMap.block(4 * t, 0, 4, 5) = this->NRGA_getT()->getTransformation();
		//outfile << this->NRGA_getT()->getTransformation();
		//outfile << std::endl;
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransformationMapFS = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3 * p_numPoints_T, 4);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> frequency = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(p_numPoints_T, 1);
	
	
	// IMPORTANT NOTE:: One doesnot know from howmany NETS or LATTICE_BLOCKS it belongs to.
	// CALCULATION of FREQUENCY HISTOGRAM  of every template point
	for (int t = 0; t < p_numPoints_T; t++)
	{
		std::vector<size_t> indexesNNT(KDTree_T->k);
		std::vector<T> distNNT(KDTree_T->k);
		KDTree_T->getKNNS(template_pts_prev[t], indexesNNT, distNNT);

		for (size_t i = 0; i < indexesNNT.size(); i++)
		{
		    TransformationMapFS.block(3 * indexesNNT.at(i), 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		    frequency(indexesNNT.at(i)) += 1;
		}
		
		TransformationMapFS.block(3 * t, 0, 3, 4) +=  TransformationMap.block(3*t, 0, 3, 4);
		frequency(t) += 1;
	}
	
	for(int t = 0; t < p_numPoints_T; t++)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Tpt_t(3,1);
		Tpt_t << template_pts_prev[t].pos.x, template_pts_prev[t].pos.y, template_pts_prev[t].pos.z;

		// NOTE: Due to this operation we get multiple goal positions 
		// ---> the solution is either "FAST SUMMATION" or "BLENDING GOAL POSITIONS" 

		Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > g_t = 
					     (((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 0, 3, 3)) * Tpt_t).transpose()				
					      + 
					     ((1.0f/frequency(t)) * TransformationMapFS.block(3 * t, 3, 3, 1)).transpose();
		
		template_points_curr[t].v.vx +=  ((g_t(0, 0) - template_pts_prev[t].pos.x))/10.0f;
		template_points_curr[t].v.vy +=  ((g_t(0, 1) - template_pts_prev[t].pos.y))/10.0f;
		template_points_curr[t].v.vz +=  ((g_t(0, 2) - template_pts_prev[t].pos.z))/10.0f;
		template_points_curr[t].kappa_vel +=  ((g_t(0, 3) - template_pts_prev[t].kappa))/5.0f;

		template_points_curr[t].v.vx = (1.1 * template_points_curr[t].v.vx) + (g_t(0, 0) - template_pts_prev[t].pos.x)/20.0f;
		template_points_curr[t].v.vy = (1.1 * template_points_curr[t].v.vy) + (g_t(0, 1) - template_pts_prev[t].pos.y)/20.0f;
		template_points_curr[t].v.vz = (1.1 * template_points_curr[t].v.vz) + (g_t(0, 2) - template_pts_prev[t].pos.z)/20.0f;
		//template_points_curr[t].kappa_vel = (1.1 * template_points_curr[t].kappa_vel) +  (g_t(0, 3) - template_pts_prev[t].kappa) /20.0f;
		
		template_points_curr[t].pos.x = template_pts_prev[t].pos.x + template_points_curr[t].v.vx * this->param->p_timeStep;
		template_points_curr[t].pos.y = template_pts_prev[t].pos.y + template_points_curr[t].v.vy * this->param->p_timeStep;
		template_points_curr[t].pos.z = template_pts_prev[t].pos.z + template_points_curr[t].v.vz * this->param->p_timeStep;
		//template_points_curr[t].kappa = template_pts_prev[t].kappa + template_points_curr[t].kappa_vel * this->param->p_timeStep;
		
		template_points_curr[t].v.vx *= 0.1;
		template_points_curr[t].v.vy *= 0.1;
		template_points_curr[t].v.vz *= 0.1;
		//template_points_curr[t].kappa_vel *= 0.1;

// 		this->NRGA_getT()->updatePointVelocity3D(template_points_curr[t].v.vx, template_points_curr[t].v.vy, template_points_curr[t].v.vz, t);
// 		this->NRGA_getT()->updatePointPosition(template_points_curr[t].pos.x, template_points_curr[t].pos.y,template_points_curr[t].pos.z,t);
	}

// 	for (int t = 0; t < p_numPoints_T; t++)
// 	{
// 		template_points_curr[t].pos.x = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.x;
// 		template_points_curr[t].pos.y = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.y;
// 		template_points_curr[t].pos.z = this->NRGA_getT()->retrieveShapePoints3D()[t].pos.z;
// 
// 		template_points_curr[t].v.vx = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vx;
// 		template_points_curr[t].v.vy = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vy;
// 		template_points_curr[t].v.vz = 0.60 * this->NRGA_getT()->retrieveShapePoints3D()[t].v.vz;
// 	}

	//outfile.close();
	TransformationMap.resize(0, 0);
	TransformationMapFS.resize(0, 0);
}

// calculate The mean Deviation(mean) between the Resulting Template and Reference points
template<class T> void NRGA<T>::NRGA_calculateMeanDeviation(point3D<T>* template_points, point3D<T>* reference_points, int numPtsTR)
{
    T deviationX = 0.0;
    T deviationY = 0.0;
    T deviationZ = 0.0;

    for (int tr = 0; tr < numPtsTR; ++tr)
    {
        deviationX = template_points[tr].pos.x - reference_points[tr].pos.x;
        deviationY = template_points[tr].pos.y - reference_points[tr].pos.y;
        deviationZ = template_points[tr].pos.z - reference_points[tr].pos.z;

        this->meanDeviationTR += std::sqrt( (deviationX * deviationX) + (deviationY * deviationY) + (deviationZ * deviationZ));
    }

    this->meanDeviationTR /= numPtsTR;
}

// calculate The Root mean Square error(RMSE) between the Resulting Template and Reference points
template<class T> void NRGA<T>::NRGA_calculateRMSE(point3D<T>* template_points, point3D<T>* reference_points, int numPtsTR, T meanDiv)
{
    T div_i      = 0.0;
    T MSE        = 0.0;
    T RMSE       = 0.0;
    T deviationX = 0.0;
    T deviationY = 0.0;
    T deviationZ = 0.0;

    for (int tr = 0; tr < numPtsTR; ++tr)
    {
        deviationX = template_points[tr].pos.x - reference_points[tr].pos.x;
        deviationY = template_points[tr].pos.y - reference_points[tr].pos.y;
        deviationZ = template_points[tr].pos.z - reference_points[tr].pos.z;

        div_i = std::sqrt( (deviationX * deviationX) + (deviationY * deviationY) + (deviationZ * deviationZ));

        MSE += (div_i - meanDiv) * (div_i - meanDiv);
    }

    MSE /= numPtsTR;

    this->RMSE_TR = std::sqrt(MSE);
}

template<class T> void NRGA<T>::NRGA_save_as_flo(std::string save_path_flo, point3D<T>* template_points_prev, point3D<T>* template_points_curr, int N)
{

    int width = 128, height = 55;

    FILE *stream = fopen(save_path_flo.c_str(), "wb");
    if (stream == 0){
        //throw CError("WriteFlowFile: could not open %s", filename);
        cout << "error while writing the file" << endl;
        exit(0);
    };

    fprintf(stream, "PIEH");

    if ((int)fwrite(&width,  sizeof(int), 1, stream) != 1 ||
    	(int)fwrite(&height, sizeof(int), 1, stream) != 1){
            cout << "error while writing the file" << endl;
            exit(0);
        };

    int nBands = 2;
    int n = nBands * width;

    for (int y1 = 0; y1 < height; y1++)
        for (int y2 = 0; y2 < width; y2++){
            float ptr_u = (float)(template_points_curr[y1*width + y2].pos.x - template_points_prev[y1*width + y2].pos.x);
            float ptr_v = (float)(template_points_curr[y1*width + y2].pos.y - template_points_prev[y1*width + y2].pos.y);
            fwrite(&ptr_u, sizeof(float), 1, stream);
            fwrite(&ptr_v, sizeof(float), 1, stream);
        }
        
    fclose(stream);
}

template<class T> void NRGA<T>::NRGA_Gravitational_Force2D(point3D<T>* template_points, point3D<T>* reference_points, point3D<T>* template_pts_prev, int p_numPoints_T, int p_numPoints_R)
{}

template<class T> void NRGA<T>::NRGA_Spring_Force2D(point3D<T>* template_points, point3D<T>* template_pts_prev, int p_numPoints_T, int p_numPoints_R)
{}

template<class T> void NRGA<T>::NRGA_Spring_Force3D(point3D<T>* template_points, point3D<T>* template_pts_prev, Eigen::Matrix<T, 3, 1>& centroid_template_curr, Eigen::Matrix<T, 3, 1>& centroid_template_prev, int p_numPoints_T, int p_numPoints_R)
{}

template<class T> void NRGA<T>::NRGA_sineDistort(point3D<T>* p, int numPts)
{}

template<class T> void NRGA<T>::plotGPE_vs_Iteration(Gnuplot& gpe, std::vector<T>& gpe_values)
{
    if (gpe_values.size() == 1)
    {
        gpe << "set xlabel 'Iterations/time'\n";
        gpe << "set ylabel 'potential Energy'\n";
        gpe << "set yrange [-1800:-1500]\n";
        gpe << "set xrange [0:400]\n";
        gpe << "set border linewidth 4\n";
        /*
        gpe << "set xlabel 'Mass'\n";
        gpe << "set ylabel 'Distance'\n";
        gpe << "set zlabel 'Gravitational Potential Energy'\n";
        gpe << "set zrange [-2.001:1]\n";
        gpe << "set yrange [-6:6]\n";
        gpe << "set xrange [0:0.00006674]\n";
        */
    }

    gpe << "plot '-' binary" << gpe.binfmt(gpe_values) << "with lines notitle linewidth 5\n";
    gpe.sendBinary(gpe_values);
    gpe.flush();
}

template<class T> NRGA<T>::~NRGA()
{
}


template class NRGA<float>;
template class NRGA<double>;
