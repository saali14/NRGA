/*
 * @author : Sk Aziz ali
 * @place  : DFKI
 * @date   : 18th May,2015
*/

#include "shape_S_points.h"

/*
 * this is the default constructor
*/
template<class T> shape_point<T>::shape_point()
{
	
    this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
    this->v.vx   = (T)0.0; this->v.vy   = (T)0.0; this->v.vz   = (T)0.0;

    this->pt3D.a = this->acc;
    this->pt3D.v = this->v;

    this->pt3D.pos.x = (T)0.0;
    this->pt3D.pos.y = (T)0.0;
    this->pt3D.pos.z = (T)0.0;

    this->pt3D.m = (T)0.0;

    this->pt3D.color.R = 0;
    this->pt3D.color.G = 0;
    this->pt3D.color.B = 0;

    this->point_X_max = (T)0.0;
    this->point_X_min = (T)0.0;
    this->point_Y_max = (T)0.0;
    this->point_Y_min = (T)0.0;
    
    // clear map of point sets for safety
    this->pointSets3D.clear();
}

template<class T> shape_point<T>::shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points,  
					      velocity<T> vel, 
					      acceleration<T> accl, 
					      p_color color, 
					      T mass)
{
    this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
    this->v.vx   = (T)0.0; this->v.vy   = (T)0.0; this->v.vz   = (T)0.0;

    this->pt3D.a = this->acc;
    this->pt3D.v = this->v;

    this->pt3D.pos.x = (T)0.0;
    this->pt3D.pos.y = (T)0.0;
    this->pt3D.pos.z = (T)0.0;

    this->pt3D.color.R = 0;
    this->pt3D.color.G = 0;
    this->pt3D.color.B = 0;

    // clear map as well 
    this->pointSets3D.clear();

    for (long num_pts = 0; num_pts < shape_points.rows(); ++num_pts) {

        this->pt3D.pos.x      = shape_points(num_pts,0);
        this->pt3D.pos.y      = shape_points(num_pts,1);
        this->pt3D.pos.z      = shape_points(num_pts,2);
        this->pt3D.a.ax		= accl.ax;
        this->pt3D.a.ay		= accl.ay;
        this->pt3D.a.az		= accl.az;
        this->pt3D.v.vx		= vel.vx;
        this->pt3D.v.vy		= vel.vy;
        this->pt3D.v.vz		= vel.vz;
	this->pt3D.color.R	= color.R;
	this->pt3D.color.G	= color.G;
	this->pt3D.color.B	= color.B;
	
	this->pt3D.m			= mass;

	// update the point set map
	this->pointSets3D.insert(std::make_pair(num_pts, this->pt3D));
    }
}

template<class T> shape_point<T>::shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points,  
					      velocity<T> vel, 
					      acceleration<T> accl, 
					      p_color color, 
					      T mass,
					      int decimation_index
 					    )
{
    this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
    this->v.vx   = (T)0.0; this->v.vy   = (T)0.0; this->v.vz   = (T)0.0;

    this->pt3D.a = this->acc;
    this->pt3D.v = this->v;

    this->pt3D.pos.x = (T)0.0;
    this->pt3D.pos.y = (T)0.0;
    this->pt3D.pos.z = (T)0.0;

    this->pt3D.color.R = 0;
    this->pt3D.color.G = 0;
    this->pt3D.color.B = 0;

    // clear map as well 
    this->pointSets3D.clear();

    for (long num_pts = 0; num_pts < shape_points.rows(); num_pts + decimation_index) {

        this->pt3D.pos.x      = shape_points(num_pts,0);
        this->pt3D.pos.y      = shape_points(num_pts,1);
        this->pt3D.pos.z      = shape_points(num_pts,2);
        this->pt3D.a.ax		= accl.ax;
        this->pt3D.a.ay		= accl.ay;
        this->pt3D.a.az		= accl.az;
        this->pt3D.v.vx		= vel.vx;
        this->pt3D.v.vy		= vel.vy;
        this->pt3D.v.vz		= vel.vz;
	this->pt3D.color.R	= color.R;
	this->pt3D.color.G	= color.G;
	this->pt3D.color.B	= color.B;
	
	this->pt3D.m		= mass;

	// update the point set map
	this->pointSets3D.insert(std::make_pair(num_pts, this->pt3D));
    }
}

template<class T> shape_point<T>::shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
					      velocity<T> vel, 
					      acceleration<T> accl, 
					      Eigen::MatrixXi shape_colors, 
					      T mass)
{
	this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
	this->v.vx = (T)0.0; this->v.vy = (T)0.0; this->v.vz = (T)0.0;

	this->pt3D.a = this->acc;
	this->pt3D.v = this->v;

	this->pt3D.pos.x = (T)0.0;
	this->pt3D.pos.y = (T)0.0;
	this->pt3D.pos.z = (T)0.0;

	this->pt3D.color.R = 0;
	this->pt3D.color.G = 0;
	this->pt3D.color.B = 0;

	// clear map as well 
	this->pointSets3D.clear();

	for (long num_pts = 0; num_pts < shape_points.rows(); ++num_pts) {
	    this->pt3D.pos.x = shape_points(num_pts, 0);
	    this->pt3D.pos.y = shape_points(num_pts, 1);
	    this->pt3D.pos.z = shape_points(num_pts, 2);
	    this->pt3D.a.ax = accl.ax;
	    this->pt3D.a.ay = accl.ay;
	    this->pt3D.a.az = accl.az;
	    this->pt3D.v.vx = vel.vx;
	    this->pt3D.v.vy = vel.vy;
	    this->pt3D.v.vz = vel.vz;

	    this->pt3D.color.R = shape_colors(num_pts, 0);
	    this->pt3D.color.G = shape_colors(num_pts, 1);
	    this->pt3D.color.B = shape_colors(num_pts, 2);

	    this->pt3D.m = mass;

	    // update the point set map
	    this->pointSets3D.insert(std::make_pair(num_pts, this->pt3D));
	}
}

template<class T> shape_point<T>::shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
					      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  curvatures,
					      velocity<T> vel, 
					      acceleration<T> accl,
					      Eigen::MatrixXi color, 
					      T mass)
{
  this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
	this->v.vx = (T)0.0; this->v.vy = (T)0.0; this->v.vz = (T)0.0;

	this->pt4D.a = this->acc;
	this->pt4D.v = this->v;

	this->pt4D.pos.x = (T)0.0;
	this->pt4D.pos.y = (T)0.0;
	this->pt4D.pos.z = (T)0.0;

	this->pt4D.color.R = 0;
	this->pt4D.color.G = 0;
	this->pt4D.color.B = 0;

	// clear map as well 
	this->pointSets4D.clear();

	for (long num_pts = 0; num_pts < shape_points.rows(); ++num_pts) {
	    this->pt4D.pos.x = shape_points(num_pts, 0);
	    this->pt4D.pos.y = shape_points(num_pts, 1);
	    this->pt4D.pos.z = shape_points(num_pts, 2);
	    this->pt4D.kappa = curvatures(num_pts, 0);
	    
	    this->pt4D.a.ax = accl.ax;
	    this->pt4D.a.ay = accl.ay;
	    this->pt4D.a.az = accl.az;
	    this->pt4D.kappa_acc = (T)0.0;
	    
	    this->pt4D.v.vx = vel.vx;
	    this->pt4D.v.vy = vel.vy;
	    this->pt4D.v.vz = vel.vz;
	    this->pt4D.kappa_vel = (T)0.0;

	    this->pt4D.color.R = color(num_pts, 0);
	    this->pt4D.color.G = color(num_pts, 1);
	    this->pt4D.color.B = color(num_pts, 2);

	    this->pt4D.m = mass;

	    // update the point set map
	    this->pointSets4D.insert(std::make_pair(num_pts, this->pt4D));
	}
}


template<class T> shape_point<T>::shape_point(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_points, 
					      velocity<T> vel, 
					      acceleration<T> accl, 
					      Eigen::MatrixXi shape_colors, 
					      T mass,
					      int decimation_index
 					    )
{
	this->acc.ax = (T)0.0; this->acc.ay = (T)0.0; this->acc.az = (T)0.0;
	this->v.vx = (T)0.0; this->v.vy = (T)0.0; this->v.vz = (T)0.0;

	this->pt3D.a = this->acc;
	this->pt3D.v = this->v;

	this->pt3D.pos.x = (T)0.0;
	this->pt3D.pos.y = (T)0.0;
	this->pt3D.pos.z = (T)0.0;

	this->pt3D.color.R = 0;
	this->pt3D.color.G = 0;
	this->pt3D.color.B = 0;

	// clear map as well 
	this->pointSets3D.clear();

	for (long num_pts = 0; num_pts < shape_points.rows(); num_pts + decimation_index) {
	    this->pt3D.pos.x = shape_points(num_pts, 0);
	    this->pt3D.pos.y = shape_points(num_pts, 1);
	    this->pt3D.pos.z = shape_points(num_pts, 2);
	    this->pt3D.a.ax = accl.ax;
	    this->pt3D.a.ay = accl.ay;
	    this->pt3D.a.az = accl.az;
	    this->pt3D.v.vx = vel.vx;
	    this->pt3D.v.vy = vel.vy;
	    this->pt3D.v.vz = vel.vz;

	    this->pt3D.color.R = shape_colors(num_pts, 0);
	    this->pt3D.color.G = shape_colors(num_pts, 1);
	    this->pt3D.color.B = shape_colors(num_pts, 2);

	    this->pt3D.m = mass;

	    // update the point set map
	    this->pointSets3D.insert(std::make_pair(num_pts, this->pt3D));
	}
}



/*
 * update shape points attributes
*/
template<class T> void shape_point<T>::updateShapePoints3D(point3D<T>* updated_Pts, long num_pts)
{
	/*
	// clear initial mapping
	this->pointSets3D.clear();

	for (long t = 0; t < num_pts; t++)
	{
		this->pt3D.pos.x = updated_Pts[t].pos.x;
		this->pt3D.pos.y = updated_Pts[t].pos.y;
		this->pt3D.pos.z = updated_Pts[t].pos.z;

		this->pt3D.a.ax = updated_Pts[t].a.ax;
		this->pt3D.a.ay = updated_Pts[t].a.ay;
		this->pt3D.a.az = updated_Pts[t].a.az;

		this->pt3D.v.vx = updated_Pts[t].v.vx;
		this->pt3D.v.vy = updated_Pts[t].v.vy;
		this->pt3D.v.vz = updated_Pts[t].v.vz;

		this->pt3D.color.R = updated_Pts[t].color.R;
		this->pt3D.color.G = updated_Pts[t].color.G;
		this->pt3D.color.B = updated_Pts[t].color.B;

		this->pt3D.m = updated_Pts[t].m;
		
		// insert pair entry into the map
		this->pointSets3D.insert(std::make_pair(t, this->pt3D));
	} // commented on 13.01.2017 
	*/
	this->pts3D = updated_Pts;
}

template<class T> void shape_point<T>::updateAcceleration3D(point3D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.a.ax = updated_Pts[t].a.ax;
		this->pointSets3D.find(t)->second.a.ay = updated_Pts[t].a.ay;
		this->pointSets3D.find(t)->second.a.az = updated_Pts[t].a.az;
	}*/ // commented on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].a.ax = updated_Pts[t].a.ax;
		this->pts3D[t].a.ay = updated_Pts[t].a.ay;
		this->pts3D[t].a.az = updated_Pts[t].a.az;
	}
}

template<class T> void shape_point<T>::updateVelocity3D(point3D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.v.vx = updated_Pts[t].v.vx;
		this->pointSets3D.find(t)->second.v.vy = updated_Pts[t].v.vy;
		this->pointSets3D.find(t)->second.v.vz = updated_Pts[t].v.vz;
	}*/ // commented on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].v.vx = updated_Pts[t].v.vx;
		this->pts3D[t].v.vy = updated_Pts[t].v.vy;
		this->pts3D[t].v.vz = updated_Pts[t].v.vz;
	}
}

template<class T> void shape_point<T>::updatePosition3D(point3D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.pos.x = updated_Pts[t].pos.x;
		this->pointSets3D.find(t)->second.pos.y = updated_Pts[t].pos.y;
		this->pointSets3D.find(t)->second.pos.z = updated_Pts[t].pos.z;
	}*/ // comment on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.x = updated_Pts[t].pos.x;
		this->pts3D[t].pos.y = updated_Pts[t].pos.y;
		this->pts3D[t].pos.z = updated_Pts[t].pos.z;
	}
}

template<class T> void shape_point<T>::updatePointVelocity3D(T UpVx, T UpVy, T UpVz, int ptIdx)
{
	this->pts3D[ptIdx].v.vx = UpVx;
	this->pts3D[ptIdx].v.vy = UpVy;
	this->pts3D[ptIdx].v.vz = UpVz;
}

template<class T> void shape_point<T>::updatePointPosition3D(T UpPx, T UpPy, T UpPz, int ptIdx)
{
	this->pts3D[ptIdx].pos.x = UpPx;
	this->pts3D[ptIdx].pos.y = UpPy;
	this->pts3D[ptIdx].pos.z = UpPz;
}

template<class T> void shape_point<T>::updateShapePoints4D(point4D<T>* updated_Pts, long num_pts)
{
	this->pts4D = updated_Pts;
}

template<class T> void shape_point<T>::updateAcceleration4D(point4D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.a.ax = updated_Pts[t].a.ax;
		this->pointSets3D.find(t)->second.a.ay = updated_Pts[t].a.ay;
		this->pointSets3D.find(t)->second.a.az = updated_Pts[t].a.az;
	}*/ // commented on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts4D[t].a.ax = updated_Pts[t].a.ax;
		this->pts4D[t].a.ay = updated_Pts[t].a.ay;
		this->pts4D[t].a.az = updated_Pts[t].a.az;
		this->pts4D[t].kappa_acc = updated_Pts[t].kappa_acc;
	}
}

template<class T> void shape_point<T>::updateVelocity4D(point4D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.v.vx = updated_Pts[t].v.vx;
		this->pointSets3D.find(t)->second.v.vy = updated_Pts[t].v.vy;
		this->pointSets3D.find(t)->second.v.vz = updated_Pts[t].v.vz;
	}*/ // commented on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts4D[t].v.vx = updated_Pts[t].v.vx;
		this->pts4D[t].v.vy = updated_Pts[t].v.vy;
		this->pts4D[t].v.vz = updated_Pts[t].v.vz;
		this->pts4D[t].kappa_vel = updated_Pts[t].kappa_vel;
	}
}

template<class T> void shape_point<T>::updatePosition4D(point4D<T>* updated_Pts)
{
	/*
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pointSets3D.find(t)->second.pos.x = updated_Pts[t].pos.x;
		this->pointSets3D.find(t)->second.pos.y = updated_Pts[t].pos.y;
		this->pointSets3D.find(t)->second.pos.z = updated_Pts[t].pos.z;
	}*/ // comment on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts4D[t].pos.x = updated_Pts[t].pos.x;
		this->pts4D[t].pos.y = updated_Pts[t].pos.y;
		this->pts4D[t].pos.z = updated_Pts[t].pos.z;
		this->pts4D[t].kappa = updated_Pts[t].kappa;
	}
}

template<class T> void shape_point<T>::updatePointVelocity4D(T UpVx, T UpVy, T UpVz, T UpVk, int ptIdx)
{
	this->pts4D[ptIdx].v.vx 	= UpVx;
	this->pts4D[ptIdx].v.vy 	= UpVy;
	this->pts4D[ptIdx].v.vz 	= UpVz;
	this->pts4D[ptIdx].kappa_vel 	= UpVk;
}

template<class T> void shape_point<T>::updatePointPosition4D(T UpPx, T UpPy, T UpPz, T UpPk, int ptIdx)
{
	this->pts4D[ptIdx].pos.x = UpPx;
	this->pts4D[ptIdx].pos.y = UpPy;
	this->pts4D[ptIdx].pos.z = UpPz;
	this->pts4D[ptIdx].kappa = UpPk;
}


template<class T> point3D<T>* shape_point<T>::retrieveShapePoints3D()
{
	return this->pts3D;
}

template<class T> point4D<T>* shape_point<T>::retrieveShapePoints4D()
{
	return this->pts4D;
}

template<class T> void shape_point<T>::linearTranslation_X(T x)
{
	/*
	for (typename std::map<long, point<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
		(*it).second.pos.x += x;*/ // commented on 13.01.2017
	
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.x += x;
	}
}

template<class T> void shape_point<T>::linearTranslation_Y(T y)
{
	/*
	for (typename std::map<long, point<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
		(*it).second.pos.y += y;*/ // commented on 13.01.2017
	
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.y += y;
	}
}

template<class T> void shape_point<T>::linearTranslation_Z(T z)
{
	/*
	for (typename std::map<long, point<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
		(*it).second.pos.z += z;*/ // commented on 13.01.2017
	
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.z += z;
	}
}

template<class T> void shape_point<T>::linearTranslation_XYZ(T x,T y, T z)
{
	/*for (typename std::map<long, point<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		(*it).second.pos.x += x;
		(*it).second.pos.y += y;
		(*it).second.pos.z += z;
	}*/ // commented on 13.01.2017

	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.x += x;
		this->pts3D[t].pos.y += y;
		this->pts3D[t].pos.z += z;
	}
}

template<class T> void shape_point<T>::translatePoint(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> translation, size_t ptIdx)
{
	this->pts3D[ptIdx].pos.x += translation(0);
	this->pts3D[ptIdx].pos.y += translation(1);
	this->pts3D[ptIdx].pos.z += translation(2);
}

template<class T> void shape_point<T>::addTranslation_XYZ(Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>& translatation)
{
	if (translatation.rows() == this->pointSets3D.size())
	{
		/*for (typename std::map<long, point<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
		{
			(*it).second.pos.x += translatation(std::distance(this->pointSets3D.begin(), it), 0);
			(*it).second.pos.y += translatation(std::distance(this->pointSets3D.begin(), it), 1);
			(*it).second.pos.z += translatation(std::distance(this->pointSets3D.begin(), it), 2);
		}*/ // commented on 13.01.2017

		for (long t = 0; t < this->pointSets3D.size(); t++)
		{
			this->pts3D[t].pos.x += translatation(t, 0);
			this->pts3D[t].pos.y += translatation(t, 1);
			this->pts3D[t].pos.z += translatation(t, 2);
		}
	}
	
}

template<class T> void shape_point<T>::rotateAbout_X_axis(T theta_x)
{
    Eigen::Matrix<T, 3, 1> Unit_X_axis;
    Unit_X_axis << 1,
                   0,
                   0;

    Eigen::Matrix<T, 3, 3> m_rotation  = Eigen::AngleAxis<T>(theta_x , Unit_X_axis).matrix();
    Eigen::Matrix<T, 3, 1> rpt;

	for (typename std::map<long, point3D<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		rpt <<	(*it).second.pos.x,
				(*it).second.pos.y,
				(*it).second.pos.z;

		rpt = m_rotation * rpt;

		(*it).second.pos.x = rpt(0, 0);
		(*it).second.pos.y = rpt(1, 0);
		(*it).second.pos.z = rpt(2, 0);
	}
}

template<class T> void shape_point<T>::rotateAbout_Y_axis(T theta_y)
{
    Eigen::Matrix<T, 3, 1> Unit_Y_axis;
    Unit_Y_axis << (T)0.0,
                   (T)1.0,
                   (T)0.0;

    Eigen::Matrix<T,3,3> m_rotation  = Eigen::AngleAxis<T>(theta_y, Unit_Y_axis).matrix();
    Eigen::Matrix<T, 3, 1> rpt;

	for (typename std::map<long, point3D<T>>::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		rpt <<	(*it).second.pos.x,
				(*it).second.pos.y,
				(*it).second.pos.z;

		rpt = m_rotation * rpt;

		(*it).second.pos.x = rpt(0, 0);
		(*it).second.pos.y = rpt(1, 0);
		(*it).second.pos.z = rpt(2, 0);
	}
}

template<class T> void shape_point<T>::rotateAbout_Z_axis(T theta_z)
{
    Eigen::Matrix<T, 3, 1> Unit_Z_axis;
    Unit_Z_axis << 0,
                   0,
                   1;

    Eigen::Matrix<T, 3, 3> m_rotation  = Eigen::AngleAxis<T>(theta_z , Unit_Z_axis).matrix();
    Eigen::Matrix<T, 3, 1> rpt;

	for (typename std::map<long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		rpt <<	(*it).second.pos.x,
				(*it).second.pos.y,
				(*it).second.pos.z;

		rpt = m_rotation * rpt;

		(*it).second.pos.x = rpt(0, 0);
		(*it).second.pos.y = rpt(1, 0);
		(*it).second.pos.z = rpt(2, 0);
	}
}

template<class T> void shape_point<T>::rotateAbout_XYZ_axis(Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> rotation)
{	
	Eigen::Matrix<T, 3, 1> rpt;

	for (typename std::map<long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		rpt <<	(*it).second.pos.x,
			(*it).second.pos.y,
			(*it).second.pos.z;

		rpt = rotation * rpt;

		(*it).second.pos.x = rpt(0, 0);
		(*it).second.pos.y = rpt(1, 0);
		(*it).second.pos.z = rpt(2, 0);
	} // commented on 13.01.2017

// 	Eigen::Matrix<T, 3, 1> rpt;
// 
// 	for (long t = 0; t < this->pointSets3D.size(); t++)
// 	{
// 		rpt << this->pts3D[t].pos.x,
// 			   this->pts3D[t].pos.y,
// 			   this->pts3D[t].pos.z;
// 
// 		rpt = rotation * rpt;
// 
// 		this->pts3D[t].pos.x = rpt(0, 0);
// 		this->pts3D[t].pos.y = rpt(1, 0);
// 		this->pts3D[t].pos.z = rpt(2, 0);
// 	}

}

template<class T> void shape_point<T>::rotateAbout_XY_axis(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation)
{
	Eigen::Matrix<T, 2, 1> rpt;

	for (typename std::map<long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		rpt <<	(*it).second.pos.x,
			(*it).second.pos.y;

		rpt = rotation * rpt;

		(*it).second.pos.x = rpt(0, 0);
		(*it).second.pos.y = rpt(1, 0);
		(*it).second.pos.z = 1.0;
	}
}

template<class T> void shape_point<T>::rotatePointXYZ(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation, size_t ptIdx)
{
	Eigen::Matrix<T, 3, 1> rpt;

	rpt << 	this->pts3D[ptIdx].pos.x,
		this->pts3D[ptIdx].pos.y,
		this->pts3D[ptIdx].pos.z;

	rpt = rotation * rpt;
	
	this->pts3D[ptIdx].pos.x = rpt(0, 0);
	this->pts3D[ptIdx].pos.y = rpt(1, 0);
	this->pts3D[ptIdx].pos.z = rpt(2, 0);
}


template<class T> void shape_point<T>::scaleShape(T sx, T sy, T sz)
{
	for (long t = 0; t < this->pointSets3D.size(); t++)
	{
		this->pts3D[t].pos.x *= sx;
		this->pts3D[t].pos.y *= sy;
		this->pts3D[t].pos.z *= sz;
	}
}

template<class T> void shape_point<T>::scalePoint(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> scale, size_t ptIdx)
{
	this->pts3D[ptIdx].pos.x *= scale(0);
	this->pts3D[ptIdx].pos.y *= scale(1);
	this->pts3D[ptIdx].pos.z *= scale(2);
}

template<class T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  shape_point<T>::getTransformationMatrix(T x,
                                                                                                            T y,
                                                                                                            T z,
                                                                                                            T theta_x,
                                                                                                            T theta_y,
                                                                                                            T theta_z,
                                                                                                            T s)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> transformationMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,4);
	return transformationMatrix;
}

template<class T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shape_point<T>::getTransformation()
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tr(this->rotation.rows(), this->rotation.cols() + this->translation.cols());
	tr << this->rotation, this->translation;

	return tr;
}

template<class T> void shape_point<T>::setTransformationMatrix(T x, T y, T z, T theta_x, T theta_y, T theta_z, T s)
{
}

template<class T> void shape_point<T>::mapToBoundingBox(T new_X_min, T new_X_max, T new_Y_min, T new_Y_max, T X_min, T X_max, T Y_min, T Y_max)
{
	for (typename std::map< long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	{
		(*it).second.pos.x = ((new_X_max - new_X_min)*((*it).second.pos.x - X_min) / (X_max - X_min)) + new_X_min;
		(*it).second.pos.y = ((new_Y_max - new_Y_min)*((*it).second.pos.y - Y_min) / (Y_max - Y_min)) + new_Y_min;
	}

	this->point_X_min = new_X_min;
	this->point_X_max = new_X_max;
	this->point_Y_min = new_Y_min;
	this->point_Y_max = new_Y_max;
}



/*
 * Normalization of point sets individually
*/
template<class T> void shape_point<T>::normalizePoints()
{
    this->point_X_min = this->point_X_max = this->pointSets3D.find(0)->second.pos.x;
    this->point_Y_min = this->point_Y_max = this->pointSets3D.find(0)->second.pos.y;
    this->point_Z_min = this->point_Z_max = this->pointSets3D.find(0)->second.pos.z;

    for (int i = 1; i < this->pointSets3D.size(); ++i) {

        // get maximum and minimum of X values
        if(this->pointSets3D.find(i)->second.pos.x > this->point_X_max)
            this->point_X_max = this->pointSets3D.find(i)->second.pos.x;
        else if(this->pointSets3D.find(i)->second.pos.x < this->point_X_min)
            this->point_X_min = this->pointSets3D.find(i)->second.pos.x;

        // get maximum and minimum of Y Values
        if(this->pointSets3D.find(i)->second.pos.y > this->point_Y_max)
            this->point_Y_max = this->pointSets3D.find(i)->second.pos.y;
        else if(this->pointSets3D.find(i)->second.pos.y < this->point_Y_min)
            this->point_Y_min = this->pointSets3D.find(i)->second.pos.y;

        // get maximum and minimum of Z Values
        if(this->pointSets3D.find(i)->second.pos.z > this->point_Z_max)
            this->point_Z_max = this->pointSets3D.find(i)->second.pos.z;
        else if(this->pointSets3D.find(i)->second.pos.z < this->point_Z_min)
            this->point_Z_min = this->pointSets3D.find(i)->second.pos.z;
    }

    std::cout << " Maximum value of X : " << this->point_X_max << " Minimum value of X : " << this->point_X_min << '\n' <<
                 " Maximum value of Y : " << this->point_Y_max << " Minimum value of Y : " << this->point_Y_min << '\n' <<
                 " Maximum value of Z : " << this->point_Z_max << " Minimum value of Z : " << this->point_Z_min << '\n' << std::endl;

    for (typename std::map<long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
    {
	(*it).second.pos.x = ((*it).second.pos.x) / (this->point_X_max - this->point_X_min);
	(*it).second.pos.y = ((*it).second.pos.y) / (this->point_Y_max - this->point_Y_min);
	if(this->point_Z_max!=this->point_Z_min)
	  (*it).second.pos.z = ((*it).second.pos.z) / (this->point_Z_max - this->point_Z_min);
	else 
	  (*it).second.pos.z = 0.0;
    }
}

template<class T> void shape_point<T>::applyMassNormalization()
{
   this->point_Mass_Max = this->point_Mass_Min = this->pointSets3D.find(0)->second.m;
   for (int i = 1; i < this->pointSets3D.size(); ++i) 
   {
      // get maximum and minimum of X values
        if(this->pointSets3D.find(i)->second.m > this->point_Mass_Max)
            this->point_Mass_Max = this->pointSets3D.find(i)->second.m;
        else if(this->pointSets3D.find(i)->second.m < this->point_Mass_Min)
            this->point_Mass_Min = this->pointSets3D.find(i)->second.m;
   }
   
    for (typename std::map<long, point3D<T> >::iterator it = this->pointSets3D.begin(); it != this->pointSets3D.end(); it++)
	(*it).second.m = ((*it).second.m) / (this->point_Mass_Max - this->point_Mass_Min);
}


/*
 * calculate the number of Noisy points as per mentioned percentage
*/
template<class T> long shape_point<T>::calcNumberofNoisyPoints(T percentage, long numPts_T)
{
    long total_point = (long)std::floor((numPts_T)/(1.0 - (percentage/100.0)));
    long noisyPts    = total_point - numPts_T;

    //std::cout << "Number of noisy points for " << percentage << "% of noise : " << noisyPts << std::endl;
    return noisyPts;
}

/*
 * Add Gaussian noise to 3D Template points
*/
template<class T> void shape_point<T>::add3DGaussianNoise(T standard_deviation, long num_ptsT, velocity<T> vel, acceleration<T> accl, p_color color, T mass)
{
    position<T> center_of_Mass;

    for (long i = 0; i < this->pointSets3D.size(); ++i)
    {
        // calculate conter of mass of the set of points
        center_of_Mass.x += this->pointSets3D.find(i)->second.pos.x;
        center_of_Mass.y += this->pointSets3D.find(i)->second.pos.y;
        center_of_Mass.z += this->pointSets3D.find(i)->second.pos.z;
    }

    center_of_Mass.x /= this->pointSets3D.size();
    center_of_Mass.y /= this->pointSets3D.size();
    center_of_Mass.z /= this->pointSets3D.size();

    // Gaussian Noise generation
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    std::normal_distribution<T> gaussian_NoiseX(center_of_Mass.x, standard_deviation);
    std::normal_distribution<T> gaussian_NoiseY(center_of_Mass.y, standard_deviation);
    std::normal_distribution<T> gaussian_NoiseZ(center_of_Mass.z, standard_deviation);


    long total_pts = this->pointSets3D.size();

    for (long i = 0; i < num_ptsT; ++i) {

        this->pt3D.pos.x = gaussian_NoiseX(generator);
        this->pt3D.pos.y = gaussian_NoiseY(generator);
        this->pt3D.pos.z = gaussian_NoiseZ(generator);

        this->pt3D.a.ax = accl.ax;
        this->pt3D.a.ay = accl.ay;
        this->pt3D.a.az = accl.az;

        this->pt3D.v.vx = vel.vx;
        this->pt3D.v.vy = vel.vy;
        this->pt3D.v.vz = vel.vz;

        this->pt3D.color.R = color.R;
        this->pt3D.color.G = color.G;
        this->pt3D.color.B = color.B;
        this->pt3D.color.A = color.A;

        this->pt3D.m = mass;

        this->pointSets3D.insert(std::make_pair(total_pts + i, this->pt3D));
    }
}

/*
 * Add Gaussian Noise to 2D Template Points
*/
template<class T> void shape_point<T>::add2DGaussianNoise(T standard_deviation, long num_ptsNoise, velocity<T> vel, acceleration<T> accl, p_color color, T mass)
{
    position<T> center_of_Mass;

    for (long i = 0; i < this->pointSets3D.size(); ++i)
    {
        // calculate conter of mass of the set of points
        center_of_Mass.x += this->pointSets3D.find(i)->second.pos.x;
        center_of_Mass.y += this->pointSets3D.find(i)->second.pos.y;
    }

    center_of_Mass.x /= this->pointSets3D.size();
    center_of_Mass.y /= this->pointSets3D.size();

    // Gaussian Noise generation
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    std::normal_distribution<T> gaussian_NoiseX(center_of_Mass.x, standard_deviation);
    std::normal_distribution<T> gaussian_NoiseY(center_of_Mass.y, standard_deviation);

    long total_pts = this->pointSets3D.size();

    for (long i = 0; i < num_ptsNoise; ++i) {

        this->pt3D.pos.x = gaussian_NoiseX(generator);
        this->pt3D.pos.y = gaussian_NoiseY(generator);
        this->pt3D.pos.z = this->pointSets3D.find(i)->second.pos.z;

        this->pt3D.a.ax = accl.ax;
        this->pt3D.a.ay = accl.ay;
        this->pt3D.a.az = accl.az;

        this->pt3D.v.vx = vel.vx;
        this->pt3D.v.vy = vel.vy;
        this->pt3D.v.vz = vel.vz;

        this->pt3D.color.R = color.R;
        this->pt3D.color.G = color.G;
        this->pt3D.color.B = color.B;
        this->pt3D.color.A = color.A;

        this->pt3D.m = mass;

        this->pointSets3D.insert(std::make_pair(total_pts + i, this->pt3D));
    }
}


/*
 * add 3D Uniformly distributed noise on 3D Template Points
*/
template<class T> void shape_point<T>::addUniformNoise3D(long num_ptsT, velocity<T> vel, acceleration<T> accl, p_color color, T mass)
{
    this->point_X_min = this->point_X_max = this->pointSets3D.find(0)->second.pos.x;
    this->point_Y_min = this->point_Y_max = this->pointSets3D.find(0)->second.pos.y;
    this->point_Z_min = this->point_Z_max = this->pointSets3D.find(0)->second.pos.z;

    for (long i = 1; i < this->pointSets3D.size(); ++i) {
        // get maximum and minimum of X values
        if(this->pointSets3D.find(i)->second.pos.x > this->point_X_max)
            this->point_X_max = this->pointSets3D.find(i)->second.pos.x;
        else if(this->pointSets3D.find(i)->second.pos.x < this->point_X_min)
            this->point_X_min = this->pointSets3D.find(i)->second.pos.x;

        // get maximum and minimum of Y Values
        if(this->pointSets3D.find(i)->second.pos.y > this->point_Y_max)
            this->point_Y_max = this->pointSets3D.find(i)->second.pos.y;
        else if(this->pointSets3D.find(i)->second.pos.y < this->point_Y_min)
            this->point_Y_min = this->pointSets3D.find(i)->second.pos.y;

        // get maximum and minimum of Z Values
        if(this->pointSets3D.find(i)->second.pos.z > this->point_Z_max)
            this->point_Z_max = this->pointSets3D.find(i)->second.pos.z;
        else if(this->pointSets3D.find(i)->second.pos.z < this->point_Z_min)
            this->point_Z_min = this->pointSets3D.find(i)->second.pos.z;
    }

    T maxXDist = this->point_X_max - this->point_X_min;
    T maxYDist = this->point_Y_max - this->point_Y_min;
    T maxZDist = this->point_Z_max - this->point_Z_min;

    T diagDist = std::sqrt((maxXDist*maxXDist) + (maxYDist*maxYDist) + (maxZDist*maxZDist));

    this->point_X_min = this->point_X_min - (std::fabs(maxXDist)* 0.2);
    this->point_Y_min = this->point_Y_min - (std::fabs(maxYDist)* 0.2);
    this->point_Z_min = this->point_Z_min - (std::fabs(maxZDist)* 0.2);
    this->point_X_max = this->point_X_max + (std::fabs(maxXDist)* 0.2);
    this->point_Y_max = this->point_Y_max + (std::fabs(maxYDist)* 0.2);
    this->point_Z_max = this->point_Z_max + (std::fabs(maxZDist)* 0.2);


    std::cout << " Maximum value of X : " << this->point_X_max << " Minimum value of X : " << this->point_X_min << '\n' <<
                 " Maximum value of Y : " << this->point_Y_max << " Minimum value of Y : " << this->point_Y_min << '\n' <<
                 " Maximum value of Z : " << this->point_Z_max << " Minimum value of Z : " << this->point_Z_min << '\n' << std::endl;


    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<T> uniformDistX(this->point_X_min , this->point_X_max);
    std::uniform_real_distribution<T> uniformDistY(this->point_Y_min , this->point_Y_max);
    std::uniform_real_distribution<T> uniformDistZ(this->point_Z_min , this->point_Z_max);

    long total_pts = this->pointSets3D.size();

    for (long i = 0; i < num_ptsT; ++i) {

        this->pt3D.pos.x = uniformDistX(generator);
        this->pt3D.pos.y = uniformDistY(generator);
        this->pt3D.pos.z = uniformDistZ(generator);

        this->pt3D.a.ax = accl.ax;
        this->pt3D.a.ay = accl.ay;
        this->pt3D.a.az = accl.az;

        this->pt3D.v.vx = vel.vx;
        this->pt3D.v.vy = vel.vy;
        this->pt3D.v.vz = vel.vz;

        this->pt3D.color.R = color.R;
        this->pt3D.color.G = color.G;
        this->pt3D.color.B = color.B;
        this->pt3D.color.A = color.A;

        this->pt3D.m = mass;

        this->pointSets3D.insert(std::make_pair(total_pts + i, this->pt3D));
    }
}

/*
 * add 2D Uniformly distributed noise on 2D Template points
*/
template<class T> void shape_point<T>::addUniformNoise2D(long num_ptsT, velocity<T> vel, acceleration<T> accl, p_color color, T mass)
{
    this->point_X_min = this->point_X_max = this->pointSets3D.find(0)->second.pos.x;
    this->point_Y_min = this->point_Y_max = this->pointSets3D.find(0)->second.pos.y;
    this->point_Z_min = this->point_Z_max = this->pointSets3D.find(0)->second.pos.z;

    for (long i = 1; i < this->pointSets3D.size(); ++i) {

        // get maximum and minimum of X values
        if(this->pointSets3D.find(i)->second.pos.x > this->point_X_max)
            this->point_X_max = this->pointSets3D.find(i)->second.pos.x;
        else if(this->pointSets3D.find(i)->second.pos.x < this->point_X_min)
            this->point_X_min = this->pointSets3D.find(i)->second.pos.x;

        // get maximum and minimum of Y Values
        if(this->pointSets3D.find(i)->second.pos.y > this->point_Y_max)
            this->point_Y_max = this->pointSets3D.find(i)->second.pos.y;
        else if(this->pointSets3D.find(i)->second.pos.y < this->point_Y_min)
            this->point_Y_min = this->pointSets3D.find(i)->second.pos.y;

        // get maximum and minimum of Z Values
        if(this->pointSets3D.find(i)->second.pos.z > this->point_Z_max)
            this->point_Z_max = this->pointSets3D.find(i)->second.pos.z;
        else if(this->pointSets3D.find(i)->second.pos.z < this->point_Z_min)
            this->point_Z_min = this->pointSets3D.find(i)->second.pos.z;
    }

    T maxXDist = this->point_X_max - this->point_X_min;
    T maxYDist = this->point_Y_max - this->point_Y_min;
    T maxZDist = this->point_Z_max - this->point_Z_min;

    T diagDist = std::sqrt((maxXDist*maxXDist) + (maxYDist*maxYDist) + (maxZDist*maxZDist));

    this->point_X_min = this->point_X_min - (std::fabs(maxXDist)* 0.2);
    this->point_Y_min = this->point_Y_min - (std::fabs(maxYDist)* 0.2);
    this->point_X_max = this->point_X_max + (std::fabs(maxXDist)* 0.2);
    this->point_Y_max = this->point_Y_max + (std::fabs(maxYDist)* 0.2);


    std::cout << " Maximum value of X : " << this->point_X_max << " Minimum value of X : " << this->point_X_min << '\n' <<
                 " Maximum value of Y : " << this->point_Y_max << " Minimum value of Y : " << this->point_Y_min << '\n' <<
                 " Maximum value of Z : " << this->point_Z_max << " Minimum value of Z : " << this->point_Z_min << '\n' << std::endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<T> uniformDistX(this->point_X_min , this->point_X_max);
    std::uniform_real_distribution<T> uniformDistY(this->point_Y_min , this->point_Y_max);

    long total_pts = this->pointSets3D.size();


    for (long i = 0; i < num_ptsT; ++i)
    {
        this->pt3D.pos.x = uniformDistX(generator);
        this->pt3D.pos.y = uniformDistY(generator);
        this->pt3D.pos.z = this->point_Z_min;

        this->pt3D.a.ax = accl.ax;
        this->pt3D.a.ay = accl.ay;
        this->pt3D.a.az = accl.az;

        this->pt3D.v.vx = vel.vx;
        this->pt3D.v.vy = vel.vy;
        this->pt3D.v.vz = vel.vz;

        this->pt3D.color.R = color.R;
        this->pt3D.color.G = color.G;
        this->pt3D.color.B = color.B;
        this->pt3D.color.A = color.A;

        this->pt3D.m = mass;

        this->pointSets3D.insert(std::make_pair(total_pts + i, this->pt3D));
    }
}

template<typename T> void shape_point<T>::applyMissingChunk(long percentage)
{  
      // We have to store the labels or indices of those missing points 
      // and inform reference points to reorder themselves by pushing them into the indices
  
    long nupPts          = this->pointSets3D.size();
    long missing_points  = (long)std::floor(nupPts * (percentage/100.0));
    long UpperboundIndex = nupPts - missing_points - 1;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<T> uniformRand(1, (T)UpperboundIndex);
    
    long startIndex = (long)std::floor(uniformRand(generator)) - 1;
    long endIndex   = startIndex +  missing_points - 1;
     
    this->missingStartIdx = startIndex; 
    this->missingEndIdx = endIndex;
    
    for (long it =  startIndex; it < endIndex; it++)
	this->pointSets3D.erase(it);
    
    for (long it = endIndex; it < nupPts; it++)
    {
         this->pt3D = this->pointSets3D.find(it)->second;
	 this->pointSets3D.erase(it);
	 this->pointSets3D.insert(std::make_pair(startIndex++, this->pt3D));
    }
}

template<typename T> void shape_point<T>::reOrderMissingIndexing(long startIdx, long endIdx)
{
    // First  --> store the missing index points in a temporary map 
    // Second --> Delete or erase the desired desired chunk from original map 
    long startIndex =  startIdx;
    long endIndex   = endIdx;
    long nupPts = this->pointSets3D.size();
    
    std::map<long, point3D<T> > TempPointSets3D;
    for (long it =  startIdx; it < endIdx; it++)
    {
	  TempPointSets3D.insert(std::make_pair(it, this->pointSets3D.find(it)->second));
	  this->pointSets3D.erase(it);
    }
    
    // restore the missing points with new Indexing
    for (long it = endIndex; it < nupPts; it++)
    {
	  this->pt3D = this->pointSets3D.find(it)->second;
	  this->pointSets3D.erase(it);
	  this->pointSets3D.insert(std::make_pair(startIndex++, this->pt3D));
    }
  
    // Add the deleted 
   for (long it = 0; it < TempPointSets3D.size(); it++)
   {
	  this->pt3D = TempPointSets3D.find(startIdx + it)->second;
	  this->pointSets3D.insert(std::make_pair(endIdx + it, this->pt3D));
   } 
}

template<typename T> long shape_point<T>::getMissingStartIdx()
{
  return this->missingStartIdx;
}

template<typename T> long shape_point<T>::getMissingEndIdx()
{
  return this->missingEndIdx;
}

template<typename T> void shape_point<T>::addPerturbations(T max_deviation)
{
    this->point_X_min = this->point_X_max = this->pointSets3D.find(0)->second.pos.x;
    this->point_Y_min = this->point_Y_max = this->pointSets3D.find(0)->second.pos.y;
    this->point_Z_min = this->point_Z_max = this->pointSets3D.find(0)->second.pos.z;

    for (long i = 1; i < this->pointSets3D.size(); ++i) {
        // get maximum and minimum of X values
        if(this->pointSets3D.find(i)->second.pos.x > this->point_X_max)
            this->point_X_max = this->pointSets3D.find(i)->second.pos.x;
        else if(this->pointSets3D.find(i)->second.pos.x < this->point_X_min)
            this->point_X_min = this->pointSets3D.find(i)->second.pos.x;

        // get maximum and minimum of Y Values
        if(this->pointSets3D.find(i)->second.pos.y > this->point_Y_max)
            this->point_Y_max = this->pointSets3D.find(i)->second.pos.y;
        else if(this->pointSets3D.find(i)->second.pos.y < this->point_Y_min)
            this->point_Y_min = this->pointSets3D.find(i)->second.pos.y;

        // get maximum and minimum of Z Values
        if(this->pointSets3D.find(i)->second.pos.z > this->point_Z_max)
            this->point_Z_max = this->pointSets3D.find(i)->second.pos.z;
        else if(this->pointSets3D.find(i)->second.pos.z < this->point_Z_min)
            this->point_Z_min = this->pointSets3D.find(i)->second.pos.z;
    }

    T maxXDist = this->point_X_max - this->point_X_min;
    T maxYDist = this->point_Y_max - this->point_Y_min;
    T maxZDist = this->point_Z_max - this->point_Z_min;

    T m_dist = std::sqrt((maxXDist*maxXDist) + (maxYDist*maxYDist) + (maxZDist*maxZDist));
    T largest_magnitude = max_deviation * m_dist;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<T> uniformRand(0, 1);
    
    for (long i = 1; i < this->pointSets3D.size(); ++i) {
      
      this->pointSets3D.find(i)->second.pos.x += uniformRand(generator) * largest_magnitude - (largest_magnitude/2 );
      this->pointSets3D.find(i)->second.pos.y += uniformRand(generator) * largest_magnitude - (largest_magnitude/2 );
      this->pointSets3D.find(i)->second.pos.z += uniformRand(generator) * largest_magnitude - (largest_magnitude/2 );
    }
}

	
template<typename T> void shape_point<T>::absOrientation(
							Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& T0,  
							Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& T1,
							size_t numPtsT0, 
							size_t numPtsT1,
							size_t idxRootPt,
							std::vector<size_t>& idxLeafPts 
							)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 1);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 1);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
	Eigen::Matrix<T, Eigen::Dynamic, 1> S;
	Eigen::Matrix<T, 3, 1> centroidT0;
	Eigen::Matrix<T, 3, 1> centroidT1;
	
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ak(3,3);
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rk(3,3);
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> Sk(3);
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tk;
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> S;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT0;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> centroidT1;

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

	//TODO: This part is for PSO related cognitive intelligence ..........................................................................

	//T0_transformed = (Rk * T0.transpose()).colwise() + Tk;			// ----> Pbest of every particle is the column vector 
	//Centroid_T0_transformed = T0_transformed.colwise().mean();			// ----> Gbest of the swarm
	//(T0_transformed.colwise() - (T0.transpose()));				// ----> (Pbest - Current) position of every particle - cognitive awareness
	//((-T0).transpose().colwise() + Centroid_T0_transformed);			// ----> (Gbest - Current) position of every particle - social awareness
	//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> r1 = 
	//	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Random(numPtsT0, 3);	// ----> random numbers for cognitive components
	//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> r2 =
	//	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Random(numPtsT0, 3);	// ----> random numbers for social components
	//....................................................................................................................................

	//for (size_t i = 0; i < idxLeafPts.size(); i++)
	//{
	//	this->rotatePointXYZ(Rk, idxLeafPts.at(i));
	//	//this->scalePoint(Sk, idxLeafPts.at(i));
	//	this->translatePoint(Tk, idxLeafPts.at(i));
	//}
	
	//this->rotatePointXYZ(Rk, idxRootPt);
	//this->scalePoint(Sk, idxRootPt);
	//this->translatePoint(Tk, idxRootPt);
	this->rotation = Rk;
	this->translation = Tk;
	//this->scale = Sk;
	

	//Ak.resize(0, 0);
	//Sk.resize(0, 0);
	//Rk.resize(0, 0);
	//Tk.resize(0, 0);
	//VT.resize(0, 0);
	//W.resize(0, 0);
	//WVT.resize(0, 0);
}

template<class T> void shape_point<T>::applyKabsch3D(point3D<T>* template_pts_in, long numPoint_Tin, point3D<T>* template_pts_out, long numPoint_Tout)
{
	T scale;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 3);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> template_mat_in;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> template_mat_out;
	Eigen::Matrix<T, Eigen::Dynamic, 1> S;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(3, 3), W(3, 3);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
	Eigen::Matrix<T, 3, 1> translation;
	Eigen::Matrix<T, 3, 1> centroid_T_in;
	Eigen::Matrix<T, 3, 1> centroid_T_out;
	template_mat_in = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPoint_Tin, 3);
	template_mat_out = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPoint_Tout, 3);

	// -----> calculation of scaling factor for template points
	T dist_template_in = 0.0;
	T dist_template_out = 0.0;

	for (long t = 0; t < numPoint_Tin - 1; t++)
	{
		dist_template_in += std::pow((template_pts_in[t + 1].pos.x - template_pts_in[t].pos.x), 2)
			+ std::pow((template_pts_in[t + 1].pos.y - template_pts_in[t].pos.y), 2)
			+ std::pow((template_pts_in[t + 1].pos.z - template_pts_in[t].pos.z), 2);
	}

	for (long r = 0; r < numPoint_Tout - 1; r++)
	{
		dist_template_out += std::pow((template_pts_out[r + 1].pos.x - template_pts_out[r].pos.x), 2)
			+ std::pow((template_pts_out[r + 1].pos.y - template_pts_out[r].pos.y), 2)
			+ std::pow((template_pts_out[r + 1].pos.z - template_pts_out[r].pos.z), 2);
	}

	if (dist_template_in <= 0.0 || dist_template_out <= 0.0)
	{
		rotation = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(3, 3);
		scale = 1.0;
		translation << 0.0, 0.0, 0.0;
	}
	else
	{
		scale = std::sqrt(dist_template_out/dist_template_in);
		template_mat_in *= scale;

		// -----> calculation of centroid of two point sets -- template points and reference points
		position<T> centroid_template_in;  centroid_template_in.x = (T)0.0; centroid_template_in.y = (T)0.0; centroid_template_in.z = (T)0.0;
		position<T> centroid_template_out; centroid_template_out.x = (T)0.0; centroid_template_out.y = (T)0.0; centroid_template_out.z = (T)0.0;

		for (long t_in = 0; t_in < numPoint_Tin; t_in++)
		{
			centroid_template_in.x += template_pts_in[t_in].pos.x;
			centroid_template_in.y += template_pts_in[t_in].pos.y;
			centroid_template_in.z += template_pts_in[t_in].pos.z;
		}

		centroid_template_in.x = centroid_template_in.x / numPoint_Tin;
		centroid_template_in.y = centroid_template_in.y / numPoint_Tin;
		centroid_template_in.z = centroid_template_in.z / numPoint_Tin;

		for (long t_out = 0; t_out < numPoint_Tout; t_out++)
		{
			centroid_template_out.x += template_pts_out[t_out].pos.x;
			centroid_template_out.y += template_pts_out[t_out].pos.y;
			centroid_template_out.z += template_pts_out[t_out].pos.z;
		}

		centroid_template_out.x = centroid_template_out.x / numPoint_Tout;
		centroid_template_out.y = centroid_template_out.y / numPoint_Tout;
		centroid_template_out.z = centroid_template_out.z / numPoint_Tout;

		// -----> shifting all the point sets so that their origin lies in the centre of mass of the respective pointsets

		for (long t_in = 0; t_in < numPoint_Tin; t_in++)
		{
			template_mat_in(t_in, 0) = template_pts_in[t_in].pos.x - centroid_template_in.x;
			template_mat_in(t_in, 1) = template_pts_in[t_in].pos.y - centroid_template_in.y;
			template_mat_in(t_in, 2) = template_pts_in[t_in].pos.z - centroid_template_in.z;
		}

		for (long t_out = 0; t_out < numPoint_Tout; t_out++)
		{
			template_mat_out(t_out, 0) = template_pts_out[t_out].pos.x - centroid_template_out.x;
			template_mat_out(t_out, 1) = template_pts_out[t_out].pos.y - centroid_template_out.y;
			template_mat_out(t_out, 2) = template_pts_out[t_out].pos.z - centroid_template_out.z;
		}

		// -----> third step is calculation of " Rotation " using SVD of covariance matrix:  transpose(template_pts) * reference_pts
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(template_mat_in.transpose() * template_mat_out, Eigen::ComputeThinU | Eigen::ComputeThinV);
		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();

		WVT = W * VT;
		S(0) = 1.0f;
		S(1) = 1.0f;
		S(2) = (WVT.determinant() <= 0) ? -1.0f : 1.0f; // correction for reflection scenario
		//std::cout << "Singular Values : " << S << std::endl;
		rotation = W * S.asDiagonal() * VT;
		//std::cout << "Rotation matrix values : "  << rotation << std::endl;

		// ------> calculation of " Translation " parameter Tx, Ty , Tz
		centroid_T_in << centroid_template_in.x,
			centroid_template_in.y,
			centroid_template_in.z;

		centroid_T_out << centroid_template_out.x,
			centroid_template_out.y,
			centroid_template_out.z;
		translation = (centroid_T_out - scale * rotation * centroid_T_in);
	}

	this->rotation = rotation;
	this->translation = translation;

	//this->scaleShape(scale, scale, scale);
	this->rotateAbout_XYZ_axis(rotation);
	//this->linearTranslation_XYZ(translation(0, 0), translation(1, 0), translation(2, 0));

	// --------> free the memory of matrices as well
	template_mat_in.resize(0, 0);
	template_mat_out.resize(0, 0);
	S.resize(0);
	VT.resize(0, 0);
	W.resize(0, 0);
	WVT.resize(0, 0);
}

template<class T> void shape_point<T>::applyKabsch2D(point3D<T>* template_pts_in, long numPoint_Tin, point3D<T>* template_pts_out, long numPoint_Tout)
{
	T scale;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rotation = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(2, 2);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> template_mat_in;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> template_mat_out;
	Eigen::Matrix<T, Eigen::Dynamic, 1> S;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> VT(2, 2), W(2, 2);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> WVT;
	Eigen::Matrix<T, 2, 1> translation; 
	Eigen::Matrix<T, 2, 1> centroid_T_in; 
	Eigen::Matrix<T, 2, 1> centroid_T_out;
	template_mat_in = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPoint_Tin, 2);
	template_mat_out = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPoint_Tout, 2);

	// -----> calculation of scaling factor for template points 
	T dist_template_in = 0.0;
	T dist_template_out = 0.0;

	for (long t = 0; t < numPoint_Tin - 1; t++)
	{
		dist_template_in += std::sqrt(std::pow((template_pts_in[t + 1].pos.x - template_pts_in[t].pos.x), 2)
			+ std::pow((template_pts_in[t + 1].pos.y - template_pts_in[t].pos.y), 2)
			);
	}

	for (long r = 0; r < numPoint_Tout - 1; r++)
	{
		dist_template_out += std::sqrt(std::pow((template_pts_out[r + 1].pos.x - template_pts_out[r].pos.x), 2)
			+ std::pow((template_pts_out[r + 1].pos.y - template_pts_out[r].pos.y), 2)
			);
	}
	
	

	if (dist_template_in <= 0.0 || dist_template_out <= 0.0)
	{
		rotation = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(2, 2);
		scale = 1.0;
		translation << 0.0, 0.0, 0.0;
	}
	else
	{
		scale = 1.0; //dist_template_out/dist_template_in ;
		template_mat_in *= scale;

		// -----> calculation of centroid of two point sets -- template points and reference points

		position<T> centroid_template_in;  centroid_template_in.x = (T)0.0; centroid_template_in.y = (T)0.0; centroid_template_in.z = (T)0.0;
		position<T> centroid_template_out; centroid_template_out.x = (T)0.0; centroid_template_out.y = (T)0.0; centroid_template_out.z = (T)0.0;

		for (long t_in = 0; t_in < numPoint_Tin; t_in++)
		{
			centroid_template_in.x += template_pts_in[t_in].pos.x;
			centroid_template_in.y += template_pts_in[t_in].pos.y;
		}

		centroid_template_in.x = centroid_template_in.x / numPoint_Tin;
		centroid_template_in.y = centroid_template_in.y / numPoint_Tin;

		for (long t_out = 0; t_out < numPoint_Tout; t_out++)
		{
			centroid_template_out.x += template_pts_out[t_out].pos.x;
			centroid_template_out.y += template_pts_out[t_out].pos.y;
		}

		centroid_template_out.x = centroid_template_out.x / numPoint_Tout;
		centroid_template_out.y = centroid_template_out.y / numPoint_Tout;

		// -----> shifting all the point sets so that their origin lies in the centre of mass of the respective pointsets

		for (long t_in = 0; t_in < numPoint_Tin; t_in++)
		{
			template_mat_in(t_in, 0) = template_pts_in[t_in].pos.x - centroid_template_in.x;
			template_mat_in(t_in, 1) = template_pts_in[t_in].pos.y - centroid_template_in.y;
		}

		for (long t_out = 0; t_out < numPoint_Tout; t_out++)
		{
			template_mat_out(t_out, 0) = template_pts_out[t_out].pos.x - centroid_template_out.x;
			template_mat_out(t_out, 1) = template_pts_out[t_out].pos.y - centroid_template_out.y;
		}

		// -----> third step is calculation of SVD of covariance matrix:  transpose(template_pts) * reference_pts

		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(template_mat_in.transpose() * template_mat_out, Eigen::ComputeThinU | Eigen::ComputeThinV);

		S = svd.singularValues();
		VT = svd.matrixU().transpose();
		W = svd.matrixV();


		//std::cout << "U values : "  << U << std::endl;
		//std::cout << "VT values : " << VT << std::endl;

		WVT = W * VT;
		S(0) = 1.0f;
		S(1) = (WVT.determinant() <= 0) ? -1.0f : 1.0f; // correction for reflection scenario

		//std::cout << "Singular Values : " << S << std::endl;
		rotation = W * S.asDiagonal() * VT;
		//std::cout << "Rotation matrix values : "  << rotation << std::endl;


		// ------> calculation of translation parameter Tx, Ty , Tz
		

		centroid_T_in << centroid_template_in.x,
			centroid_template_in.y;

		centroid_T_out << centroid_template_out.x,
			centroid_template_out.y;

		translation = (centroid_T_out - scale * rotation * centroid_T_in);
		
	}
	
	scaleShape(scale, scale, 1.0);
	rotateAbout_XY_axis(rotation);
	linearTranslation_XYZ(translation(0, 0), translation(1, 0), 0.0);
	

	// free the memory of matrices as well
	template_mat_in.resize(0, 0);
	template_mat_out.resize(0, 0);
	S.resize(0);
	VT.resize(0, 0);
	W.resize(0, 0);
	WVT.resize(0, 0);
}

template<class T> shape_point<T>::~shape_point()
{
	// map is cleared
	this->pointSets3D.clear();
}

template class shape_point<float> ;
template class shape_point<double> ;
