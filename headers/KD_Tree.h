/* @author : Sk Aziz Ali
*  @place  : German Research Center for Artificial Intelligence, Germany 
*  @date   : 14th June, 2015
*/


#ifndef _KDTREE_ 
#define _KDTREE_ 

#pragma once
// C++ Standard
#include <ctime>					// clock
#include <cmath>					// math routines
#include <cstring>					// C string ops
#include <fstream>					// file I/O
#include <memory>					// Memory 
#include <cstdlib>

// Approximate Nearest Neighbourhood 
#include <nanoflann.hpp>

// Program headers
#include "shape_S_points.h"

using namespace nanoflann;
using namespace std;
using namespace Eigen;

#ifndef CLOCKS_PER_SEC					// define clocks-per-second if needed
  #define CLOCKS_PER_SEC 1000000
#endif

// Structure builds KD tree on "Normal Splitting" and Serches "Nearest Neighbourhood" for a query point using " NANO-FLANN library "
template <typename T>
struct kdt_nanoFlann
{
	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >  my_kd_tree_t;

	int nPts;				// total number of points 
	const size_t  k;			// number of nearest neighbors
	size_t  dim = 3;			// dimension
	T eps;					// error bound
	int maxPts;				// maximum number of data points
	point3D<T>* pts;			// points from Direct Memory Access chunk out of on/off chip GPU  	
	my_kd_tree_t* kdtree;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat;

	kdt_nanoFlann(point3D<T>* points, int numPts, size_t dimensionality, T errorPrecision, size_t numNeighbours) : 
		pts(points), nPts(numPts), dim(dimensionality), eps(errorPrecision), k(numNeighbours)
	{
		// map the points only to an Eigen matrix 
		//this->mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(nPts, dim);
		this->mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				&points[0].pos.x,
				nPts,
				dim,
				Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point<T>, m) / sizeof(T) + 1, 1));

		// build the kd tree and their nodal indices
		this->kdtree = new my_kd_tree_t(this->dim, this->mat, 10 /* max leaf */);
		this->kdtree->index->buildIndex();
	}

	void getKNNS(point3D<T> qPoint, std::vector<size_t> &ret_indexes, std::vector<T> &out_dists_sqr)
	{
		// map the query point to  stl vector
		std::vector<T> qpt(dim);
		qpt[0] = qPoint.pos.x; qpt[1] = qPoint.pos.y; qpt[2] = qPoint.pos.z;

		// do a knn search
		nanoflann::KNNResultSet<T> resultSet(this->k);
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		this->kdtree->index->findNeighbors(resultSet, &qpt[0], nanoflann::SearchParams(10));
	}

	void printNeighbours(std::vector<size_t> &ret_indexes, std::vector<T> &out_dists_sqr)
	{
		std::cout << "knnSearch(nn=" << k << "): \n";
		for (size_t i = 0; i<k; i++)
			std::cout << "ret_index[" << i << "]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;
	}

	~kdt_nanoFlann()
	{
		delete[] this->kdtree;
	}
};



template <typename T>
struct kdt_nanoFlann4D
{
	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >  my_kd_tree_t;

	int nPts;				// total number of points 
	const size_t  k;			// number of nearest neighbors
	size_t  dim = 4;			// dimension
	T eps;					// error bound
	int maxPts;				// maximum number of data points
	point4D<T>* pts;			// points from Direct Memory Access chunk out of on/off chip GPU  	
	my_kd_tree_t* kdtree;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat;

	kdt_nanoFlann4D(point4D<T>* points, int numPts, size_t dimensionality, T errorPrecision, size_t numNeighbours) : 
		pts(points), nPts(numPts), dim(dimensionality), eps(errorPrecision), k(numNeighbours)
	{
		// map the points only to an Eigen matrix 
		//this->mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(nPts, dim);
		this->mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				&points[0].pos.x,
				nPts,
				dim,
				Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<T>, m) / sizeof(T) + 1, 1));

		// build the kd tree and their nodal indices
		this->kdtree = new my_kd_tree_t(this->dim, this->mat, 10 /* max leaf */);
		this->kdtree->index->buildIndex();
	}

	void getKNNS(point4D<T> qPoint, std::vector<size_t> &ret_indexes, std::vector<T> &out_dists_sqr)
	{
		// map the query point to  stl vector
		std::vector<T> qpt(dim);
		qpt[0] = qPoint.pos.x; qpt[1] = qPoint.pos.y; qpt[2] = qPoint.pos.z; qpt[0] = qPoint.kappa;

		// do a knn search
		nanoflann::KNNResultSet<T> resultSet(this->k);
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		this->kdtree->index->findNeighbors(resultSet, &qpt[0], nanoflann::SearchParams(10));
	}

	void printNeighbours(std::vector<size_t> &ret_indexes, std::vector<T> &out_dists_sqr)
	{
		std::cout << "knnSearch(nn=" << k << "): \n";
		for (size_t i = 0; i<k; i++)
			std::cout << "ret_index[" << i << "]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;
	}

	~kdt_nanoFlann4D()
	{
		delete[] this->kdtree;
	}
};
#endif
