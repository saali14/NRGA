/*
* @author: Sk Aziz Ali
* @place : DFKI
*/

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>
#include <tuple>
#include <math.h>
#include <map>
#include <utility>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

struct p_color
{
	unsigned char R; unsigned char G; unsigned char B; unsigned char A;
};

template<typename T>
struct position
{
	T x; T y; T z;
};

template<typename T>
struct velocity
{
	T vx; T vy; T vz;
};

template<typename T>
struct acceleration
{
	T ax; T ay; T az;
};

template<typename T>
struct point
{
	position<T>	pos;    // position of the shape point
	velocity<T>	v;      // velocity of shape point
	acceleration<T>	a;	// acceleration of shape point
	p_color		color;	// color of particle or shape points
	T		m;	// mass of particles/points
};

template<typename T>
struct point3D
{
	position<T>	pos;	// Position 
	velocity<T>	v;	// Velocity 
	acceleration<T>	a;      // Acceleration 
	p_color		color;	// color of particle or shape points
	T		m;	// Mass
};

template<typename T>
struct point4D
{
	position<T>	pos;	// Position 3D
	T 		kappa;	// Curvature
	
	velocity<T>	v;      // Velocity 
	T		kappa_vel;
	
	acceleration<T>	a;      // Acceleration
	T		kappa_acc;
	
	p_color		color;	// color of particle or shape points
	T		m;	// Mass
};

template<typename T>
struct spring{

};

template<typename T>
class shape_point
{
public:
	T                                                   m;
	velocity<T>                                         v;
	acceleration<T>                                     acc;
	p_color                                             color;
	point3D<T>                                          pt3D;
	point3D<T>*                                         pts3D;
	point4D<T>                                          pt4D;
	point4D<T>*                                         pts4D;
	std::map<long, point3D<T> >                         pointSets3D;
	std::map<long, point4D<T> >                         pointSets4D;
	T                                                   point_X_min;
	T                                                   point_X_max;
	T                                                   point_Y_min;
	T                                                   point_Y_max;
	T                                                   point_Z_min;
	T                                                   point_Z_max;
	T						    point_Mass_Max;
	T						    point_Mass_Min;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    rotation;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    scale;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    translation;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    all_rotations;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    all_scales;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>    all_translation;
	long missingStartIdx, missingEndIdx;
	
	shape_point();
	shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
		    velocity<T> vel, 
		    acceleration<T> accl, 
		    p_color color, 
		    T mass);
	shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
		    velocity<T> vel, 
		    acceleration<T> accl, 
		    p_color color, 
		    T mass, 
		    int decimation_index
 		  );
	shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
		    velocity<T> vel, 
		    acceleration<T> accl, 
		    Eigen::MatrixXi color, 
		    T mass);
	shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
		    velocity<T> vel, 
		    acceleration<T> accl, 
		    Eigen::MatrixXi color, 
		    T mass,
		    int decimation_index
 		  );
	shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
		    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  curvatures,
		    velocity<T> vel, 
		    acceleration<T> accl,
		    Eigen::MatrixXi color, 
		    T mass);

	point3D<T>*						retrieveShapePoints3D();
	point4D<T>*						retrieveShapePoints4D();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllRotations();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllTranslations();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllScales();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllRotations4D();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllTranslations4D();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	retrieveAllScales4D();
	
	void updateShapePoints3D(point3D<T>*, long);
	void updateVelocity3D(point3D<T>*);
	void updateAcceleration3D(point3D<T>*);
	void updatePosition3D(point3D<T>*);
	void updateColor3D(point3D<T>*);
	void updateMass3D(point3D<T>*);
	void updatePointVelocity3D(T UpVx, T UpVy, T UpVz, int ptIdx);
	void updatePointPosition3D(T UpPx, T UpPy, T UpPz, int ptIdx);
	
	void updateShapePoints4D(point4D<T>*, long);
	void updateVelocity4D(point4D<T>*);
	void updateAcceleration4D(point4D<T>*);
	void updatePosition4D(point4D<T>*);
	void updateColor4D(point4D<T>*);
	void updateMass4D(point4D<T>*);
	void updatePointVelocity4D(T UpVx, T UpVy, T UpVz, T UpVk, int ptIdx);
	void updatePointPosition4D(T UpPx, T UpPy, T UpPz, T UpPk, int ptIdx);

	void linearTranslation_X(T x);
	void linearTranslation_Y(T y);
	void linearTranslation_Z(T z);
	void linearTranslation_XYZ(T x, T y, T z);
	void translatePoint(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> translation, size_t ptIdx);
	void addTranslation_XYZ(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& translatation);

	void rotateAbout_X_axis(T theta_x);
	void rotateAbout_Y_axis(T theta_y);
	void rotateAbout_Z_axis(T theta_z);
	void rotateAbout_XYZ_axis(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation);
	void rotateAbout_XY_axis(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation);
	void rotatePointXYZ(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation, size_t ptIdx);
	void scaleShape(T sx, T sy, T sz);
	void scalePoint(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> scale, size_t ptIdx);
	void setTransformationMatrix(T x, T y, T z, T theta_x, T theta_y, T theta_z, T s);
	void mapToBoundingBox(T new_X_min, T new_X_max, T new_Y_min, T new_Y_max, T X_min, T X_max, T Y_min, T Y_max);
	void normalizePoints();
	void absOrientation(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& T0, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& T1, size_t numPtsTO, size_t numPtsT1, size_t idxRootPt, std::vector<size_t>& idxLeafPts);
	void applyKabsch3D(point3D<T>*, long, point3D<T>*, long);
	void applyKabsch2D(point3D<T>*, long, point3D<T>*, long);
	void addUniformNoise3D(long, velocity<T> vel, acceleration<T> accl, p_color color, T mass);
	void addUniformNoise2D(long, velocity<T> vel, acceleration<T> accl, p_color color, T mass);
	void add2DGaussianNoise(T standard_deviation, long, velocity<T> vel, acceleration<T> accl, p_color color, T mass);
	void add3DGaussianNoise(T standard_deviation, long, velocity<T> vel, acceleration<T> accl, p_color color, T mass);
	void addPerturbations(T max_deviation);
	void applyMissingChunk(long);
	void applySubSampling3D();
	void applyMassNormalization();
	void calcScaleFactor(point<T>*, long, point<T>*, long);
	void build_KDTree(point<T>*, long);
	void reOrderMissingIndexing(long startIdx, long endIdx);
	long calcNumberofNoisyPoints(T percentage, long numPts_T);
	

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getCurrentShapeMatrix();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getTransformedPositionMatrix(T sx, T sy, T sz, T tx, T ty, T tz, Eigen::Matrix<T, 3, 3> rotation);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getTransformationMatrix(T x, T y, T z, T theta_x, T theta_y, T theta_z, T s);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getTransformedPoints(point<T>* initPoint, T x, T y, T z, T theta_x, T theta_y, T theta_z, T scale, velocity<T> v, acceleration<T> a);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getTransformation();
        long getMissingStartIdx();
	long getMissingEndIdx();
	virtual ~shape_point();
};
