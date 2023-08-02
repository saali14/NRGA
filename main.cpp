
/*
 * @Author : Sk Aziz Ali
 * @place  : DFKI GmbH, Kaiserslautern, Germany
 * 
 */

// Qt
#define GLEW_STATIC
#include <QApplication>
#include <qnamespace.h>
#include <qmenu.h>
#include <qcursor.h>
#include <qimage.h>
#include <qdatetime.h>
#include <GL/glew.h>
#include <QtOpenGL/QGLWidget>
#include <QDebug>
#include <QEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QTimer>
#include <QSurfaceFormat>
#include <QThread>
#include <QPainter>
#include <QColor>
#include <QFontMetrics>
#include <QRect>
#include <QBrush>
#include <QFont>
#include <QImage>
#include <QPen>

// C++11
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <ctime>
#include <time.h>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <iterator>
#include <random>
#include <algorithm>
#include <chrono>
#include <utility>

#ifdef __linux__
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#ifdef _WIN32
#include <stddef.h>
#endif

#ifdef __linux__
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#endif

// CUDA
//#include <cuda.h>
//#include <cuda_runtime.h>

// PROGRAM
//#include "cuda/add_vect.cuh"
#include "GA.h"
#include "NRGA.h"
#include "shaders.h"
#include "ply_Module_ali.h"

// BOOST
#include <boost/lambda/lambda.hpp>

void test();
float testTimeElapsed(int N);
void appendStatistics(std::string path, std::string mean , std::string RMSE);


using namespace std;

class GLWidget : public QGLWidget
{
    QTimer* timer = new QTimer(this);
    
public:
    //Parameterized constructor of GlWidget Class 
    GLubyte*            	Imagedata;
    NRGA<double>*  		NRGA_sim = NULL;
    GA<double>* 		GA_sim = NULL;
    shaders*       		NRGA_shadersR = NULL;
    shaders*   		NRGA_shadersT = NULL;
    ply_module*		NRGA_plyModuleT = NULL;
    ply_module*    		NRGA_plyModuleR = NULL;
    point3D<double>* 	template_pointsRef3D = NULL;	// template points 3D as reference 	
    point3D<double>* 	Tpts3D = NULL; 	                // template points for Direct Memory Access(DMA) in CPU
    point3D<double>* 	Rpts3D = NULL; 	                // reference points for Direct Memory Access(DMA) in GPU
    point3D<double>* 	Tpts3D_prev = NULL; 	        // Initial template points for Direct Memory Access(DMA) in CPU
    point3D<double>* 	template_pts_GPU3D = NULL; 	// template points for Direct Memory Access(DMA) in CPU
    point3D<double>* 	reference_pts_GPU3D = NULL; 	// reference points for Direct Memory Access(DMA) in GPU
    point4D<double>* 	template_points4D = NULL; 	// template points for Direct Memory Access(DMA) in CPU
    point4D<double>* 	reference_points4D = NULL; 	// reference points for Direct Memory Access(DMA) in GPU
    point4D<double>* 	template_pts_prev4D = NULL; 	// Initial template points for Direct Memory Access(DMA) in CPU
    point4D<double>* 	template_pts_GPU4D = NULL; 	// template points for Direct Memory Access(DMA) in CPU
    point4D<double>* 	reference_pts_GPU4D = NULL; 	// reference points for Direct Memory Access(DMA) in GPU
    int*                	IndicesT = NULL;
    int*                	IndicesR = NULL;
    double*             	normalsT = NULL;
    double*             	normalsR = NULL;
    char*               	colOffseT = NULL;
    char*               	colOffseR = NULL;
    char*               	accOffseT = NULL;
    char*               	velOffseT = NULL;
    char*               	massOffseT = NULL;
    char*			curvatureOffsetT = NULL;
    char*			curvatureOffsetR = NULL;
    
    kdt_nanoFlann4D<double>* KDTree4D_T = NULL;
    kdt_nanoFlann4D<double>* KDTree4D_R = NULL;
    kdt_nanoFlann<double>* 	KDTree3D_T = NULL;
    kdt_nanoFlann<double>* 	KDTree3D_R = NULL;
    
    QString infoTxt =   QString("Key G : Toggle drawing of grid") + QString("\n") +
    QString("Key I : Take snapshot") + QString("\n") +
    QString("Key S : Toggle Simulation Start/Stop") + QString("\n") +
    QString("Key P : Save Ply of Template Object") + QString("\n") +
    QString("Key H : See Help information");
    
    GLWidget(NRGA<double>* nrga, 
             GA<double>* ga, 
             shaders* shaderR, 
             shaders* shaderT, 
             ply_module* plyModuleT, 
             ply_module* plyModuleR) :
             NRGA_sim(nrga),
             GA_sim(ga),
             NRGA_shadersR(shaderR),
             NRGA_shadersT(shaderT),
             NRGA_plyModuleT(plyModuleT),
             NRGA_plyModuleR(plyModuleR)
             {
                 // shape points color
                 p_color colorT = { (unsigned char)0, (unsigned char)0, (unsigned char)255, (unsigned char)1 };
                 p_color colorR = { (unsigned char)255, (unsigned char)0, (unsigned char)0, (unsigned char)1 };
                 
                 velocity<double> vel;
                 vel.vx = 0.0f;
                 vel.vy = 0.0f;
                 vel.vz = 0.0f;
                 
                 acceleration<double>  accl;
                 accl.ax = 0.0f;
                 accl.ay = 0.0f;
                 accl.az = 0.0f;
                 
                 double p_mass = 1.0f;
                 
                 shape_point<double> template_pts;
                 shape_point<double> reference_pts;
                 
                 if(this->NRGA_sim->NRGA_getParameters()->dimT == 3)
                 {
                     if(NRGA_plyModuleT->getColors().size() == 0)
                         template_pts = shape_point<double>(NRGA_plyModuleT->getVertices(), vel, accl, colorT, p_mass);
                     else
                         template_pts = shape_point<double>(NRGA_plyModuleT->getVertices(), vel, accl, NRGA_plyModuleT->getColors(), p_mass);
                 }
                 else if(this->NRGA_sim->NRGA_getParameters()->dimT == 4)
                 {
                     if(NRGA_plyModuleT->getColors().size() > 0 && NRGA_plyModuleT->getCurvatures().size() > 0)
                         template_pts = shape_point<double>(NRGA_plyModuleT->getVertices(), NRGA_plyModuleT->getCurvatures(), vel, accl, NRGA_plyModuleT->getColors(), p_mass);
                     else
                         std::cout << "Vertex Attributes is NOT initialized for TEMPLATE" << std::endl;
                 }
                 if(this->NRGA_sim->NRGA_getParameters()->p_normalize_flag == 1)
                 {
                     template_pts.normalizePoints();
                 }
                 
                 if(this->NRGA_sim->NRGA_getParameters()->p_mass_normalize_flag == 1)
                     template_pts.applyMassNormalization();
                 
                 
                 // To apply Missing chunk 
                 //template_pts.applyMissingChunk(60);
                 
                 if(this->NRGA_sim->NRGA_getParameters()->p_noise_types == "G" 
                     && this->NRGA_sim->NRGA_getParameters()->dimT == 3)
                 {
                     
                     template_pts.add2DGaussianNoise(this->NRGA_sim->NRGA_getParameters()->p_noiseG_variance, 
                                                     template_pts.calcNumberofNoisyPoints( this->NRGA_sim->NRGA_getParameters()->p_noise_percentage,
                                                                                           template_pts.pointSets3D.size()
                                                     ), 
                                                     vel, 
                                                     accl, 
                                                     colorT, 
                                                     p_mass);
                 }
                 if(this->NRGA_sim->NRGA_getParameters()->p_noise_types == "U" 
                     && this->NRGA_sim->NRGA_getParameters()->dimT == 3)
                 {
                     
                     template_pts.addUniformNoise2D( template_pts.calcNumberofNoisyPoints( this->NRGA_sim->NRGA_getParameters()->p_noise_percentage,
                                                                                           template_pts.pointSets3D.size()
                     ), 
                     vel, 
                     accl, 
                     colorT, 
                     p_mass);
                 }
                 if(this->NRGA_sim->NRGA_getParameters()->p_noise_types == "P" 
                     && this->NRGA_sim->NRGA_getParameters()->dimT == 3)
                 {
                     double step_size = 0.05; 
                     double max_deviation = step_size;
                     template_pts.addPerturbations(max_deviation);
                 }
                 
                 // 		unsigned seedR = std::chrono::system_clock::now().time_since_epoch().count();
                 // 		std::mt19937 generatorR(seedR);
                 // 		std::uniform_real_distribution<double> uniformDistRx(M_PI/5.0 , M_PI/5.0);
                 // 		std::uniform_real_distribution<double> uniformDistRy(-M_PI/5.0 , M_PI/5.0);
                 // 		std::uniform_real_distribution<double> uniformDistRz(-M_PI/5.0 , M_PI/5.0);
                 // 		
                 // 		
                 // 		Eigen::AngleAxisd rollAngle(uniformDistRx(generatorR),  Eigen::Vector3d::UnitX());
                 // 		Eigen::AngleAxisd pitchAngle(0.0, Eigen::Vector3d::UnitY());
                 // 		Eigen::AngleAxisd yawAngle(uniformDistRz(generatorR),  Eigen::Vector3d::UnitZ());
                 // 		Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle;
                 //  		Eigen::Matrix3d rotationMatrix;
                 // 		rotationMatrix << 0.96f, -0.08f, -0.27f, 
                 // 				 -0.04f, 0.91f,  -0.42f,
                 // 				  0.27f, 0.41f,   0.87f;
                 //  		template_pts.rotateAbout_XYZ_axis(rotationMatrix);
                 
                 
                 if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
                 {
                     if(NRGA_plyModuleR->getColors().size() == 0)
                         reference_pts = shape_point<double>(NRGA_plyModuleR->getVertices(), vel, accl, colorR, p_mass);
                     else
                         reference_pts = shape_point<double>(NRGA_plyModuleR->getVertices(), vel, accl, NRGA_plyModuleR->getColors(), p_mass);
                 }
                 else if(this->NRGA_sim->NRGA_getParameters()->dimR == 4)
                 {
                     if(NRGA_plyModuleR->getColors().size() > 0 && NRGA_plyModuleR->getCurvatures().size() > 0)	
                         reference_pts = shape_point<double>(NRGA_plyModuleR->getVertices(), NRGA_plyModuleR->getCurvatures(), vel, accl, NRGA_plyModuleR->getColors(), p_mass);
                     else
                         std::cout << "Vertex Attributes is NOT initialized for REFERENCE" << std::endl;
                 }
                 if(this->NRGA_sim->NRGA_getParameters()->p_normalize_flag == 1)
                 {
                     reference_pts.normalizePoints();
                 }		
                 if(this->NRGA_sim->NRGA_getParameters()->p_mass_normalize_flag == 1)
                     reference_pts.applyMassNormalization();
                 
                 
                 // template and reference points attributes initialization starts here...
                 if(this->NRGA_sim->NRGA_getParameters()->dimT == 3)
                     this->initShapePoints_T3D(template_pts);
                 else
                     this->initShapePoints_T4D(template_pts);
                 
                 if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
                     this->initShapePoints_R3D(reference_pts);
                 else
                     this->initShapePoints_R4D(reference_pts);
                 
                 this->initFaceIndices_T(NRGA_plyModuleT->getFaces());
                 this->initFaceIndices_R(NRGA_plyModuleR->getFaces());
                 this->initVerticesNormals_T(NRGA_plyModuleT->getNormals());
                 this->initVerticesNormals_R(NRGA_plyModuleR->getNormals());
                 nrga->NRGA_setT(&template_pts);
                 ga->GA_setT(&template_pts);
                 
                 if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
                     this->KDTree3D_R = new kdt_nanoFlann<double>(this->Rpts3D, p_numPoints_R, this->NRGA_sim->NRGA_getParameters()->dimR, 0.001, this->NRGA_sim->NRGA_getParameters()->k_R);
                 else
                     this->KDTree4D_R = new kdt_nanoFlann4D<double>(this->reference_points4D, p_numPoints_R, this->NRGA_sim->NRGA_getParameters()->dimR, 0.001, this->NRGA_sim->NRGA_getParameters()->k_R);
                 
                 connect(timer, &QTimer::timeout, this, &GLWidget::on_Simulation);
                 timer->setInterval(50);
             }
             
             
    //----------------------- member variuable used------------------------------------
    long int iter = 0;
    long	 frame = 0;
    int      frame_number = 0;
    int      time_elapsed;
    int      timebase = 0;
    int      fps;   
    bool     R_Color = false;               // R pressed modifier
    bool     G_Color = false;               // G pressed modifier
    bool     B_Color = false;               // B pressed modifier
    bool     initGL = false;                // initialize openGL flag
    bool     create_VBO = false;            // create Vertex Buffer Object after opnGL initialization
    bool	 onMouse = false;		
    bool     draw_KNN = false;
    bool     draw_shape = false;            // draw shape on openGL screen
    bool	 draw_gravityField = false;	// flag to draw gravitational field
    bool     draw_Image = false;		// flag to draw shape points
    bool	 draw_grid = false;		// flag to draw a grid
    bool	 draw_axis = false;		// flag to draw axis
    bool	 draw_Simulation = false;       // flag to strat simulation
    bool	 draw_Reference = true;		// flad to toggle drawing Regerence Body
    bool	 draw_Template = true;		// flad to toggle drawing Template Body
    bool	 restart = false; 		// flag to restart the system of Registration
    bool	 save_ply = false;		// flag to save ply files
    bool	 toggleCM = false;		// flag to toggle cloud Vs mesh visualization
    bool     normalize;                     // flag to normalize point4D sets
    bool	 light_flag = true;		// flag to turn on and off light for shading 
    int      window_width = 740;		// Size of openGL window width when it starts in terms of pixels
    int      window_height = 1440;		// size of openGL window height when it starts in terms of pixels
    int      window_depth = 1000;		// Depth of OpenGL Window when it starts in terms of pixels
    long	 p_numPoints_T = 0;		// Number of shape points in Template
    long	 p_numPoints_R = 0;		// Number of shape points in Reference
    long     p_numIndices_T = 0;		// Number of indices in Index buffers for faces of Template
    long     p_numIndices_R = 0;		// Number of indices in Index buffers for faces of Reference
    double   p_scale_param;			// scaling parameter of Shape model
    double   p_theta_x;			// Angle of rotation of shape model with respect to X- axis
    double   p_theta_y;			// Angle of rotation of shape model with respect to Y- axis
    double   p_theta_z;			// Angle of rotation of shape model with respect to Z- axis
    double   p_translation_x;		// Translation parameter with respect to X-axis
    double   p_translation_y;		// Translation parameter with respect to Y-axis
    double   p_translation_z;		// Translation parameter with respect to Z-axis
    double	 p_mass;    			// particle mass - same for all particles
    
    QImage   image;
    GLint    p_left_bound;			// Left orthographich projection window bound
    GLint	 p_right_bound;			// Right orhtographic projection window bound
    GLint	 p_top_bound;			// Top orthographic projection window bound
    GLint	 p_bottom_bound;		// Bottom orthographic projection window bound
    GLint	 p_nDepth_bound;		// Top orthographic projection window bound
    GLint	 p_fDepth_bound;		// Bottom orthographic projection window bound
    GLuint   p_vaoT;			// vertex Array Object for Template
    GLuint   p_vaoR;			// vertex Array Object for Reference
    GLuint	 p_vboT;			// vertex buffer object for Template points
    GLuint	 p_vboR;			// vertex buffer object for Reference points
    GLuint	 p_iboT;			// index buffer object for Template points
    GLuint	 p_iboR;			// index buffer object for Reference points
    GLuint	 p_nboT;			// normal buffer object for Template points
    GLuint	 p_nboR;			// normal buffer object for Reference points
    GLuint   p_vboT_trans;			// transform feedback buffer object for Template points
    GLuint   p_vboR_trans;			// transform feedback buffer object for Reference points
    GLuint   p_textureId;			// texture id of colormap
    GLuint	 vao = 0;			// vertex array object
    GLuint   fbo;				// frame buffer object
    GLuint   rboColor;			// render buffer object for color
    GLuint   rboDepth;			// Bottom orthographic projection window bound
    double   epsilon = 1.0;
    double	 TRACKBALL_RADIUS = 0.6f;
    double	 radius_ = -1.01;		// --> -1.01
    QPoint   last_point_2D_;
    bool	 last_point_ok_;
    Eigen::Vector3d center_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d last_point_3D_;
    GLfloat projection_matrix_[16], modelview_matrix_[16];
    
    // Camera position
    double x = 0.0;
    double y = 0.0;
    double z = 550.0;                        // initially 5 units south of origin
    double deltaMove = 0.0;                  // initially camera doesn't move
    
    // Camera direction
    double lx = 0.0, ly = 0.0, lz = -550.0; // camera points initially along z-axis
    double upX = 0.0, upY = 1.0, upZ = 0.0; // angle of rotation for the camera direction
    double deltaAngleX = 0.0;		// additional angle change when dragging
    double deltaAngleY = 0.0; 		// additional angle change when dragging
    double deltaAngleZ = 0.0;  		// additional angle change when dragging
    double FOV = 0.09f;                    // field of view for camera 0.9
    
    // initialize template points 
    void initShapePoints_T4D(shape_point<double> &template_pt)
    {
        p_numPoints_T = template_pt.pointSets4D.size();
        
        // initial template points
        this->template_pts_prev4D = (point4D<double>*)malloc(p_numPoints_T * sizeof(point4D<double>));
        this->template_points4D = (point4D<double>*)malloc(p_numPoints_T * sizeof(point4D<double>));
        
        for (int i = 0; i < p_numPoints_T; i++)
        {
            this->template_points4D[i].pos.x = template_pt.pointSets4D.find(i)->second.pos.x;
            this->template_points4D[i].pos.y = template_pt.pointSets4D.find(i)->second.pos.y;
            this->template_points4D[i].pos.z = template_pt.pointSets4D.find(i)->second.pos.z;
            this->template_points4D[i].kappa = template_pt.pointSets4D.find(i)->second.kappa;
            
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
    void initShapePoints_R4D(shape_point<double> &reference_pt)
    {
        p_numPoints_R = reference_pt.pointSets4D.size();
        
        this->reference_points4D = (point4D<double>*)malloc(p_numPoints_R * sizeof(point4D<double>));
        
        for (int i = 0; i < p_numPoints_R; i++)
        {
            this->reference_points4D[i].pos.x = reference_pt.pointSets4D.find(i)->second.pos.x;
            this->reference_points4D[i].pos.y = reference_pt.pointSets4D.find(i)->second.pos.y;
            this->reference_points4D[i].pos.z = reference_pt.pointSets4D.find(i)->second.pos.z;
            this->reference_points4D[i].kappa = reference_pt.pointSets4D.find(i)->second.kappa;
            
            std::cout << this->reference_points4D[i].kappa << std::endl;
            
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
    
    void initShapePoints_T3D(shape_point<double> &template_pt)
    {
        p_numPoints_T = template_pt.pointSets3D.size();
        
        // initial template points
        this->template_pointsRef3D = (point3D<double>*)malloc(p_numPoints_T * sizeof(point3D<double>));
        this->Tpts3D_prev = (point3D<double>*)malloc(p_numPoints_T * sizeof(point3D<double>));
        this->Tpts3D = (point3D<double>*)malloc(p_numPoints_T * sizeof(point3D<double>));
        
        for (int i = 0; i < p_numPoints_T; i++)
        {
            this->template_pointsRef3D[i].pos.x = template_pt.pointSets3D.find(i)->second.pos.x;
            this->template_pointsRef3D[i].pos.y = template_pt.pointSets3D.find(i)->second.pos.y;
            this->template_pointsRef3D[i].pos.z = template_pt.pointSets3D.find(i)->second.pos.z;
            //this->template_pointsRef3D[i].pos.z = 0.0;
            
            this->template_pointsRef3D[i].v.vx = template_pt.pointSets3D.find(i)->second.v.vx;
            this->template_pointsRef3D[i].v.vy = template_pt.pointSets3D.find(i)->second.v.vy;
            this->template_pointsRef3D[i].v.vz = template_pt.pointSets3D.find(i)->second.v.vz;
            
            this->template_pointsRef3D[i].a.ax = template_pt.pointSets3D.find(i)->second.a.ax;
            this->template_pointsRef3D[i].a.ay = template_pt.pointSets3D.find(i)->second.a.ay;
            this->template_pointsRef3D[i].a.az = template_pt.pointSets3D.find(i)->second.a.az;
            
            this->template_pointsRef3D[i].color.R = template_pt.pointSets3D.find(i)->second.color.R;
            this->template_pointsRef3D[i].color.G = template_pt.pointSets3D.find(i)->second.color.G;
            this->template_pointsRef3D[i].color.B = template_pt.pointSets3D.find(i)->second.color.B;
            
            //this->template_pointsRef3D[i].m = template_pt.pointSets3D.find(i)->second.pos.z; //template_pt.pointSets3D.find(i)->second.m;
            this->template_pointsRef3D[i].m = template_pt.pointSets3D.find(i)->second.m;
        }
        
        for (int i = 0; i < p_numPoints_T; i++)
        {
            this->Tpts3D[i].pos.x = template_pt.pointSets3D.find(i)->second.pos.x;
            this->Tpts3D[i].pos.y = template_pt.pointSets3D.find(i)->second.pos.y;
            this->Tpts3D[i].pos.z = template_pt.pointSets3D.find(i)->second.pos.z;
            //this->Tpts3D[i].pos.z = 0.0;
            
            this->Tpts3D[i].v.vx = template_pt.pointSets3D.find(i)->second.v.vx;
            this->Tpts3D[i].v.vy = template_pt.pointSets3D.find(i)->second.v.vy;
            this->Tpts3D[i].v.vz = template_pt.pointSets3D.find(i)->second.v.vz;
            
            this->Tpts3D[i].a.ax = template_pt.pointSets3D.find(i)->second.a.ax;
            this->Tpts3D[i].a.ay = template_pt.pointSets3D.find(i)->second.a.ay;
            this->Tpts3D[i].a.az = template_pt.pointSets3D.find(i)->second.a.az;
            
            this->Tpts3D[i].color.R = template_pt.pointSets3D.find(i)->second.color.R;
            this->Tpts3D[i].color.G = template_pt.pointSets3D.find(i)->second.color.G;
            this->Tpts3D[i].color.B = template_pt.pointSets3D.find(i)->second.color.B;
            
            //this->Tpts3D[i].m = template_pt.pointSets3D.find(i)->second.pos.z; //template_pt.pointSets3D.find(i)->second.m;
            this->Tpts3D[i].m = template_pt.pointSets3D.find(i)->second.m;
        }
        
    }
    
    // initialize reference points
    void initShapePoints_R3D(shape_point<double> &reference_pt)
    {
        p_numPoints_R = reference_pt.pointSets3D.size();
        
        this->Rpts3D = (point3D<double>*)malloc(p_numPoints_R * sizeof(point3D<double>));
        
        for (int i = 0; i < p_numPoints_R; i++)
        {
            this->Rpts3D[i].pos.x = reference_pt.pointSets3D.find(i)->second.pos.x;
            this->Rpts3D[i].pos.y = reference_pt.pointSets3D.find(i)->second.pos.y;
            this->Rpts3D[i].pos.z = reference_pt.pointSets3D.find(i)->second.pos.z;
            
            this->Rpts3D[i].v.vx = reference_pt.pointSets3D.find(i)->second.v.vx;
            this->Rpts3D[i].v.vy = reference_pt.pointSets3D.find(i)->second.v.vy;
            this->Rpts3D[i].v.vz = reference_pt.pointSets3D.find(i)->second.v.vz;
            
            this->Rpts3D[i].a.ax = reference_pt.pointSets3D.find(i)->second.a.ax;
            this->Rpts3D[i].a.ay = reference_pt.pointSets3D.find(i)->second.a.ay;
            this->Rpts3D[i].a.az = reference_pt.pointSets3D.find(i)->second.a.az;
            
            this->Rpts3D[i].color.R = reference_pt.pointSets3D.find(i)->second.color.R;
            this->Rpts3D[i].color.G = reference_pt.pointSets3D.find(i)->second.color.G;
            this->Rpts3D[i].color.B = reference_pt.pointSets3D.find(i)->second.color.B;
            
            this->Rpts3D[i].m = reference_pt.pointSets3D.find(i)->second.m;			
        }
    }
    
    // Initialize indices for template points
    void initFaceIndices_T(Eigen::MatrixXi facesT)
    {	  
        if (facesT.rows() != 0 && facesT.cols() != 0)
        {
            this->IndicesT = (int*)malloc(facesT.rows() * facesT.cols() * sizeof(int));
            Eigen::Map <Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(&this->IndicesT[0], facesT.rows(), facesT.cols()) = facesT;
            this->p_numIndices_T = facesT.rows() * facesT.cols();
            
            if (((sizeof(this->IndicesT)) / sizeof((this->IndicesT)[0])) == (facesT.rows() * facesT.cols()))
                std::cout << " Mapping is succesful for TEMPLATE Index buffer " << std::endl;
            else
            {
                std::cout << " WARNING: Mapping is NOT succesful for TEMPLATE Index buffer " << std::endl;
                std::cout << " To be mapped : " << (facesT.rows() * facesT.cols()) << std::endl;
                std::cout << " Actually mapped : " << ((sizeof(this->IndicesT)) / sizeof((this->IndicesT)[0])) << std::endl;
            }
        }
    }
    
    // Initialize indices for template points
    void initFaceIndices_R(Eigen::MatrixXi facesR)
    {
        if (facesR.rows() != 0 && facesR.cols() != 0)
        {
            this->IndicesR = (int*)malloc(facesR.rows() * facesR.cols() * sizeof(int));
            Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(this->IndicesR, facesR.rows(), facesR.cols()) = facesR.cast<int>();
            
            this->p_numIndices_R = facesR.rows() * facesR.cols();
            
            if ((sizeof(this->IndicesR) / sizeof(this->IndicesR[0])) == (facesR.rows() * facesR.cols()))
                std::cout << " Mapping is succesful for Index REFERENCE buffer " << std::endl;
            else
            {
                std::cout << " WARNING: Mapping is NOT succesful for REFERENCE Index buffer " << std::endl;
                std::cout << " To be mapped : " << (facesR.rows() * facesR.cols()) << std::endl;
                std::cout << " Actually mapped : " << (sizeof(this->IndicesR) / sizeof(this->IndicesR[0])) << std::endl;
            }
        }
    }
    
    // Initialize Template vertex Normals
    void initVerticesNormals_T(Eigen::MatrixXd normalsMat_T)
    {
        if (normalsMat_T.rows() != 0 && normalsMat_T.cols() != 0)
        {
            this->normalsT = (double*)malloc(normalsMat_T.rows() * normalsMat_T.cols() * sizeof(double));
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(this->normalsT, normalsMat_T.rows(), normalsMat_T.cols()) = normalsMat_T;
            
            
            if ((sizeof(this->normalsT) / sizeof(this->normalsT[0])) == (normalsMat_T.rows() * normalsMat_T.cols()))
                std::cout << " Mapping is succesful for TEMPLATE Normal buffer " << std::endl;
            else
            {
                std::cout << " WARNING: Mapping is NOT succesful for TEMPLATE Normal buffer " << std::endl;
                std::cout << " To be mapped : " << (normalsMat_T.rows() * normalsMat_T.cols()) << std::endl;
                std::cout << " Actually mapped : " << (sizeof(this->normalsT) / sizeof(this->normalsT[0])) << std::endl;
            }
        }
    }
    
    // Initialize Reference vertex Normals
    void initVerticesNormals_R(Eigen::MatrixXd normalsMat_R)
    {
        if (normalsMat_R.rows() != 0 && normalsMat_R.cols() != 0)
        {
            this->normalsR = (double*)malloc(normalsMat_R.rows() * normalsMat_R.cols() * sizeof(double));
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(this->normalsR, normalsMat_R.rows(), normalsMat_R.cols()) = normalsMat_R;
            
            if ((sizeof(this->normalsR) / sizeof(this->normalsR[0])) == (normalsMat_R.rows() * normalsMat_R.cols()))
                std::cout << " Mapping is succesful for REFERENCE Normal buffer " << std::endl;
            else
            {
                std::cout << " WARNING: Mapping is NOT succesful for REFERENCE Normal buffer " << std::endl;
                std::cout << " To be mapped : " << (normalsMat_R.rows() * normalsMat_R.cols()) << std::endl;
                std::cout << " Actually mapped : " << (sizeof(this->normalsR) / sizeof(this->normalsR[0])) << std::endl;
            }
        }
    }
    
    void set1DTexture()
    {
        // generate the specified number of texture objects 
        glGenTextures(1, &this->p_textureId);
        glBindTexture(GL_TEXTURE_1D, p_textureId);
        
        // tells OpenGL how the data that is going to be uploaded is aligned
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        
        float data[] =
        {
            0.0f, 0.0f, 1.0f, 1.0f,
            0.0f, 1.0f, 0.0f, 1.0f,
            1.0f, 0.0f, 0.0f, 1.0f
        };
        
        glTexImage1D(
            GL_TEXTURE_1D,      // Specifies the target texture. Must be GL_TEXTURE_1D or GL_PROXY_TEXTURE_1D.
            0,                  // Specifies the level-of-detail number. Level 0 is the base image level. Level n is the nth mipmap reduction image.
            GL_RGBA32F,
            3,
            0,                  // border: This value must be 0.
            GL_RGBA,
            GL_FLOAT,
            data
        );
        
        // texture sampling/filtering operation.
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_1D, 0);
    }
    
    // Initialize for openGL context creattion
    void initializeGL(){
        
        // get context opengl-version information and print
        qDebug() << "Widget OpenGl Version			: " << format().majorVersion() << "." << format().minorVersion();
        qDebug() << "Context valid fLag				: " << context()->isValid();
        qDebug() << "Really used OpenGl				: " << context()->format().majorVersion() << "." << context()->format().minorVersion();
        qDebug() << "OpenGl information: VENDOR			: " << (const char*)glGetString(GL_VENDOR);
        qDebug() << "                    RENDERDER		: " << (const char*)glGetString(GL_RENDERER);
        qDebug() << "                    VERSION		: " << (const char*)glGetString(GL_VERSION);
        qDebug() << "                    GLSL VERSION		: " << (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);
        
        /// In modelview hand is at origin
        this->initGL = !this->initGL;
        glewExperimental = GL_TRUE;
        glClearColor(this->NRGA_sim->NRGA_getParameters()->bgColor[0], this->NRGA_sim->NRGA_getParameters()->bgColor[1], this->NRGA_sim->NRGA_getParameters()->bgColor[2], this->NRGA_sim->NRGA_getParameters()->bgColor[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearDepth(1.0f);
        glEnable(GL_DEPTH_TEST);
        
        // Fog
        GLfloat fogColor[4] = { 0.3, 0.3, 0.4, 1.0 };
        glFogi(GL_FOG_MODE, GL_LINEAR);
        glFogfv(GL_FOG_COLOR, fogColor);
        glFogf(GL_FOG_DENSITY, 0.35);
        glHint(GL_FOG_HINT, GL_DONT_CARE);
        glFogf(GL_FOG_START, 5.0f);
        glFogf(GL_FOG_END, 25.0f);
        glDisable(GL_DOUBLEBUFFER);
        glEnableClientState(GL_VERTEX_ARRAY);
        glDepthFunc(GL_LEQUAL);
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glPointSize(this->NRGA_sim->NRGA_getParameters()->pointSize);
        
        
        // scene pos and size
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glGetFloatv(GL_MODELVIEW_MATRIX, this->modelview_matrix_);
        set_scene_pos(Eigen::Vector3d(0.0f, 0.0f, 0.0f), 1.0);
        
        
        // setting up the shaders ATTRIBUTES' POSITION for Reference Object 
        if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
        {	
            /// @COMMENTS:: If Vertex indices are provided --> use "MESH" based sahders else "POINT" based shaders
            if (this->p_numIndices_R > 0)
            {
                this->NRGA_shadersR->compileShaders(
                    this->NRGA_sim->param->shaderPath_RMeshVert.c_str(),
                                                    this->NRGA_sim->param->shaderPath_RMeshFrag.c_str()
                );
                
                this->NRGA_shadersR->setNormalAttribute(glGetAttribLocation(this->NRGA_shadersR->getProgShaderId(), "normal"));
                if (this->NRGA_shadersR->getNormalAttribute() < 0) {
                    std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
                }
            }
            else
            {
                this->NRGA_shadersR->compileShaders(
                    this->NRGA_sim->param->shaderPath_RPointsVert.c_str(),
                                                    this->NRGA_sim->param->shaderPath_RPointsFrag.c_str()
                );
            }
            
            this->NRGA_shadersR->setPositionUniform(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "p"));
            if (this->NRGA_shadersR->getPositionUniform() < 0) {
                std::cerr << "Shader did not contain the 'p' uniform." << std::endl;
            }
            this->NRGA_shadersR->setPositionAttribute(glGetAttribLocation(this->NRGA_shadersR->getProgShaderId(), "position"));
            if (this->NRGA_shadersR->getPositionAttribute() < 0) {
                std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
            }
            this->NRGA_shadersR->setColorAttribute(glGetAttribLocation(this->NRGA_shadersR->getProgShaderId(), "colour"));
            if (this->NRGA_shadersR->getColorAttribute() < 0) {
                std::cerr << "Shader did not contain the 'colour' attribute." << std::endl;
            }
            
            this->NRGA_shadersR->setProjectionMatrixUniform(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "MP"));
            if (this->NRGA_shadersR->getProjectionMatrixUniform() < 0) {
                std::cerr << "Shader did not contain the 'MP - projection' attribute." << std::endl;
            }
            this->NRGA_shadersR->setModelViewMatrixUniform(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "MV"));
            if (this->NRGA_shadersR->getModelViewMatrixUniform() < 0) {
                std::cerr << "Shader did not contain the 'MV - modelview' attribute." << std::endl;
            }
        }
        
        // setting up the shaders ATTRIBUTES' POSITION for Template Object 
        if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
        {
            /// @COMMENTS:: If Vertex indices are provided --> use "MESH" based sahders else "POINT" based shaders
            if (this->p_numIndices_T > 0)
            {
                this->NRGA_shadersT->compileShaders(this->NRGA_sim->param->shaderPath_TMeshVert.c_str(),
                                                    this->NRGA_sim->param->shaderPath_TMeshFrag.c_str());
                
                this->NRGA_shadersT->setNormalAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "normal"));
                if (this->NRGA_shadersT->getNormalAttribute() < 0) {
                    std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
                }
            }
            else
            {
                this->NRGA_shadersT->compileShaders(this->NRGA_sim->param->shaderPath_TPointsVert.c_str(),
                                                    this->NRGA_sim->param->shaderPath_TPointsFrag.c_str());
            }
            
            this->NRGA_shadersT->setPositionUniform(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "p"));
            if (this->NRGA_shadersT->getPositionUniform() < 0) {
                std::cerr << "Shader did not contain the 'p' uniform." << std::endl;
            }
            this->NRGA_shadersT->setPositionAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "position"));
            if (this->NRGA_shadersT->getPositionAttribute() < 0) {
                std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
            }
            this->NRGA_shadersT->setColorAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "colour"));
            if (this->NRGA_shadersT->getColorAttribute() < 0) {
                std::cerr << "Shader did not contain the 'colour' attribute." << std::endl;
            }
            this->NRGA_shadersT->setAccelerationAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "acceleration"));
            if (this->NRGA_shadersT->getAccelerationAttribute() < 0) {
                std::cerr << "Shader did not contain the 'acceleration' attribute." << std::endl;
            }
            this->NRGA_shadersT->setVelocityAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "velocity"));
            if (this->NRGA_shadersT->getVelocityAttribute() < 0) {
                std::cerr << "Shader did not contain the 'velocity' attribute." << std::endl;
            }
            this->NRGA_shadersT->setMassAttribute(glGetAttribLocation(this->NRGA_shadersT->getProgShaderId(), "mass"));
            if (this->NRGA_shadersT->getMassAttribute() < 0) {
                std::cerr << "Shader did not contain the 'mass' attribute." << std::endl;
            }
            this->NRGA_shadersT->setProjectionMatrixUniform(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "MP"));
            if (this->NRGA_shadersT->getProjectionMatrixUniform() < 0) {
                std::cerr << "Shader did not contain the 'MP - projection' attribute." << std::endl;
            }
            this->NRGA_shadersT->setModelViewMatrixUniform(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "MV"));
            if (this->NRGA_shadersT->getModelViewMatrixUniform() < 0) {
                std::cerr << "Shader did not contain the 'MV - modelview' attribute." << std::endl;
            }
        }
    }
    
    void qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
    {
        const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
        const GLdouble ymin = -ymax;
        const GLdouble xmin = ymin * aspect;
        const GLdouble xmax = ymax * aspect;
        glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
    }
    
    /// @note camera decides renderer size
    void resizeGL(int width, int height){
        
        this->window_width = width;
        this->window_height = height;
        this->window_depth = width;
        
        if (height == 0) height = 1;
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        qgluPerspective(this->FOV, (GLfloat)width / (GLfloat)height, 0.001f, width / 10.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        this->updateGL();
    }
    
    void paintGL()
    {
        //qgluPerspective(this->FOV, (GLfloat)this->window_width / (GLfloat)this->window_height, 0.1f, this->window_width / 10.0f);
        glClearColor(this->NRGA_sim->NRGA_getParameters()->bgColor[0], this->NRGA_sim->NRGA_getParameters()->bgColor[1], this->NRGA_sim->NRGA_getParameters()->bgColor[2], this->NRGA_sim->NRGA_getParameters()->bgColor[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(this->projection_matrix_);
        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(this->modelview_matrix_);
        glPointSize(this->NRGA_sim->NRGA_getParameters()->pointSize);
        
        if (!this->initGL)
        {
            // first initialize openGL
            this->initializeGL();
            this->initGL = !this->initGL;
        }
        // create VBO
        if (!this->create_VBO)
        {
            if (glewInit() == GLEW_OK)
            {
                set1DTexture();
                createVBO();
            }
            else
            {
                std::cerr << "Failed to initialize GLEW.." << std::endl;
                std::exit(-1);
            }
            this->create_VBO = !this->create_VBO;
        }
        
        // drawing the shape 
        if (this->draw_shape)
        {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            visualize();
        }
        
        if (this->draw_grid)
        {
            drawgrid();
        }
        
        if(this->draw_gravityField)
        {
            drawForceField();
        }
        if(this->draw_KNN)
        {
            drawKNN();
        }
        if (this->draw_Simulation)
        {
            timer->start();
        }
        else
        {
            timer->stop();
        }
        
        glFlush();
    }
    
    void on_Simulation()
    {
        if(iter++ <= this->NRGA_sim->NRGA_getParameters()->p_max_iterations)
        {
            this->draw_shape = !this->draw_shape;		
            if(this->NRGA_sim->NRGA_getParameters()->dimT == 3)
            {
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboT);
                glBufferData(GL_ARRAY_BUFFER, p_numPoints_T * sizeof(point3D<double>), NULL, GL_DYNAMIC_DRAW);
                this->Tpts3D = (point3D<double>*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
                this->KDTree3D_T = new kdt_nanoFlann<double>(this->Tpts3D, p_numPoints_T, this->NRGA_sim->NRGA_getParameters()->dimT, 0.001, this->NRGA_sim->NRGA_getParameters()->k_T);
                this->NRGA_sim->NRGA_setCurrIterNo(iter);
                this->NRGA_sim->NRGA_Gravitational_Force3D_ARAP3(this->Tpts3D, this->Rpts3D, this->Tpts3D_prev, this->KDTree3D_T, this->KDTree3D_R, p_numPoints_T, p_numPoints_R);
                if (this->save_ply == true) this->saveTemplate();
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else if(this->NRGA_sim->NRGA_getParameters()->dimT == 4)
            {
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboT);
                glBufferData(GL_ARRAY_BUFFER, p_numPoints_T * sizeof(point4D<double>), NULL, GL_DYNAMIC_DRAW);
                this->template_points4D = (point4D<double>*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
                this->KDTree4D_T = new kdt_nanoFlann4D<double>(this->template_points4D, p_numPoints_T, this->NRGA_sim->NRGA_getParameters()->dimT, 0.001, this->NRGA_sim->NRGA_getParameters()->k_T);
                this->NRGA_sim->NRGA_setCurrIterNo(iter);
                this->NRGA_sim->NRGA_CurveSpaceGravitational2_Force4D(this->template_points4D, this->reference_points4D, this->template_pts_prev4D, this->KDTree4D_T, this->KDTree4D_R, p_numPoints_T, p_numPoints_R);
                if (this->save_ply == true) this->saveTemplate();
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else
            {
                ;
            }
            this->draw_shape = !this->draw_shape;
            if (this->draw_Image == true) this->slotSnapshot();
            updateGL();
        }
    }
    
    // create vertex buffer object
    void createVBO()
    {
        glGenVertexArrays(1, &this->p_vaoR);
        glBindVertexArray(this->p_vaoR);
        
        if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
        {
            if (this->Rpts3D != NULL && this->p_numPoints_R > 0)
            {
                glGenBuffers(1, &this->p_vboR);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboR);
                glBufferData(GL_ARRAY_BUFFER, this->p_numPoints_R * sizeof(point3D<double>), this->Rpts3D, GL_DYNAMIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getPositionAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getPositionAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)0);
                    std::cout << "POSITION Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (&this->Rpts3D[0].color != NULL  && p_numPoints_R > 0)
            {
                this->colOffseR = (char*)((char*)& this->Rpts3D[0].color - (char*)& this->Rpts3D[0]);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getColorAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getColorAttribute(), 3, GL_UNSIGNED_BYTE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)this->colOffseR);
                    std::cout << "COLOR Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (this->normalsR != NULL && this->p_numPoints_R > 0)
            {
                glGenBuffers(1, &this->p_nboR);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_nboR);
                glBufferData(GL_ARRAY_BUFFER, 3 * p_numPoints_R * sizeof(double), this->normalsR, GL_STATIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getNormalAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getNormalAttribute(), 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);
                    std::cout << "NORMAL Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (this->IndicesR != NULL && this->p_numIndices_R > 0)
            {
                glGenBuffers(1, &this->p_iboR);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboR);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->p_numIndices_R * sizeof(int), this->IndicesR, GL_STATIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                std::cout << "INDEX Buffer is set correctly for Reference" << std::endl;
            }
            
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
        else if(this->NRGA_sim->NRGA_getParameters()->dimR == 4)
        {
            if (this->reference_points4D != NULL && p_numPoints_R > 0)
            {
                glGenBuffers(1, &this->p_vboR);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboR);
                glBufferData(GL_ARRAY_BUFFER, p_numPoints_R * sizeof(point4D<double>), this->reference_points4D, GL_DYNAMIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getPositionAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getPositionAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)0);
                    std::cout << "POSITION Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (&this->reference_points4D[0].color != NULL  && p_numPoints_R > 0)
            {
                this->colOffseR = (char*)((char*)& this->reference_points4D[0].color - (char*)& this->reference_points4D[0]);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getColorAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getColorAttribute(), 3, GL_UNSIGNED_BYTE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)this->colOffseR);
                    std::cout << "COLOR Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (this->normalsR != NULL && this->p_numPoints_R > 0)
            {
                glGenBuffers(1, &this->p_nboR);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_nboR);
                glBufferData(GL_ARRAY_BUFFER, 3 * p_numPoints_R * sizeof(double), this->normalsR, GL_STATIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersR->getNormalAttribute());
                    glVertexAttribPointer(this->NRGA_shadersR->getNormalAttribute(), 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);
                    std::cout << "NORMAL Buffer is set correctly for Reference" << std::endl;
                }
            }
            
            if (this->IndicesR != NULL && this->p_numIndices_R > 0)
            {
                glGenBuffers(1, &this->p_iboR);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboR);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->p_numIndices_R * sizeof(int), this->IndicesR, GL_STATIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                std::cout << "INDEX Buffer is set correctly for Reference" << std::endl;
            }
            
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
        else
        {
            ;
        }
        
        glGenVertexArrays(1, &this->p_vaoT);
        glBindVertexArray(this->p_vaoT);
        
        if(this->NRGA_sim->NRGA_getParameters()->dimT == 3)
        {
            if (this->Tpts3D != NULL && this->p_numPoints_T > 0)
            {
                glGenBuffers(1, &this->p_vboT);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboT);
                glBufferData(GL_ARRAY_BUFFER, this->p_numPoints_T * sizeof(point3D<double>), this->Tpts3D, GL_DYNAMIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getPositionAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getPositionAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)0);
                    std::cout << "POSITION Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->Tpts3D[0].color != NULL  && this->p_numPoints_T > 0)
            {
                this->colOffseT = (char*)((char*)& this->Tpts3D[0].color - (char*)& Tpts3D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getColorAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getColorAttribute(), 3, GL_UNSIGNED_BYTE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)this->colOffseT);
                    std::cout << "COLOUR Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->Tpts3D[0].a != NULL  && this->p_numPoints_T > 0)
            {
                this->accOffseT = (char*)((char*)& this->Tpts3D[0].a - (char*)& Tpts3D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getAccelerationAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getAccelerationAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)this->accOffseT);
                    std::cout << "ACCELERATION Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->Tpts3D[0].v != NULL  && this->p_numPoints_T > 0)
            {
                this->velOffseT = (char*)((char*)& this->Tpts3D[0].v - (char*)& Tpts3D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getVelocityAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getVelocityAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)this->velOffseT);
                    std::cout << "VELOCITY Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->Tpts3D[0].m != NULL  && this->p_numPoints_T > 0)
            {
                this->massOffseT = (char*)((char*)& this->Tpts3D[0].color - (char*)& Tpts3D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getMassAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getMassAttribute(), 1, GL_DOUBLE, GL_FALSE, sizeof(point3D<double>), (const GLvoid *)this->massOffseT);
                    std::cout << "MASS Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (this->normalsT != NULL && this->p_numPoints_T > 0)
            {
                glGenBuffers(1, &this->p_nboT);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_nboT);
                glBufferData(GL_ARRAY_BUFFER, 3 * p_numPoints_T * sizeof(double), this->normalsT, GL_STATIC_DRAW);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getNormalAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getNormalAttribute(), 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);
                    std::cout << "NORMAL Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (this->IndicesT != NULL && this->p_numIndices_T > 0)
            {
                glGenBuffers(1, &this->p_iboT);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboT);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->p_numIndices_T * sizeof(int), this->IndicesT, GL_STATIC_DRAW);
                std::cout << "INDICES Buffer is set correctly for Template" << std::endl;
            }
        }
        else if(this->NRGA_sim->NRGA_getParameters()->dimT == 4)
        {
            if (this->template_points4D != NULL && this->p_numPoints_T > 0)
            {
                glGenBuffers(1, &this->p_vboT);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_vboT);
                glBufferData(GL_ARRAY_BUFFER, p_numPoints_T * sizeof(point4D<double>), this->template_points4D, GL_DYNAMIC_DRAW);
                
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getPositionAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getPositionAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)0);
                    std::cout << "POSITION Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->template_points4D[0].color != NULL  && this->p_numPoints_T > 0)
            {
                this->colOffseT = (char*)((char*)& this->template_points4D[0].color - (char*)& template_points4D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getColorAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getColorAttribute(), 3, GL_UNSIGNED_BYTE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)this->colOffseT);
                    std::cout << "COLOUR Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->template_points4D[0].a != NULL  && this->p_numPoints_T > 0)
            {
                this->accOffseT = (char*)((char*)& this->template_points4D[0].a - (char*)& template_points4D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getAccelerationAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getAccelerationAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)this->accOffseT);
                    std::cout << "ACCELERATION Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->template_points4D[0].v != NULL  && this->p_numPoints_T > 0)
            {
                this->velOffseT = (char*)((char*)& this->template_points4D[0].v - (char*)& template_points4D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getVelocityAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getVelocityAttribute(), 3, GL_DOUBLE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)this->velOffseT);
                    std::cout << "VELOCITY Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (&this->template_points4D[0].m != NULL  && this->p_numPoints_T > 0)
            {
                this->massOffseT = (char*)((char*)& this->template_points4D[0].color - (char*)& template_points4D[0]);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getMassAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getMassAttribute(), 1, GL_DOUBLE, GL_FALSE, sizeof(point4D<double>), (const GLvoid *)this->massOffseT);
                    std::cout << "MASS Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (this->normalsT != NULL && this->p_numPoints_T > 0)
            {
                glGenBuffers(1, &this->p_nboT);
                glBindBuffer(GL_ARRAY_BUFFER, this->p_nboT);
                glBufferData(GL_ARRAY_BUFFER, 3 * p_numPoints_T * sizeof(double), this->normalsT, GL_STATIC_DRAW);
                if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
                {
                    glEnableVertexAttribArray(this->NRGA_shadersT->getNormalAttribute());
                    glVertexAttribPointer(this->NRGA_shadersT->getNormalAttribute(), 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);
                    std::cout << "NORMAL Buffer is set correctly for Template" << std::endl;
                }
            }
            
            if (this->IndicesT != NULL && this->p_numIndices_T > 0)
            {
                glGenBuffers(1, &this->p_iboT);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboT);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->p_numIndices_T * sizeof(int), this->IndicesT, GL_STATIC_DRAW);
                std::cout << "INDICES Buffer is set correctly for Template" << std::endl;
            }
        }
        else
        {}
        
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

	// delete vertex buffer object when it is of no use
	void deleteVBO()
	{
		// deletes VBO for template and buffer points
		glDeleteBuffers(1, &this->p_vboT);
		glDeleteBuffers(1, &this->p_vboR);
		glDeleteBuffers(1, &this->p_vaoT);
		glDeleteBuffers(1, &this->p_vaoR);
		glDeleteBuffers(1, &this->p_iboT);
		glDeleteBuffers(1, &this->p_iboR);
		glDeleteBuffers(1, &this->p_nboT);
		glDeleteBuffers(1, &this->p_nboR);
	}

	// create render buffer for color
	void createRBOC(void)
	{
		glGenRenderbuffers(1, &rboColor);
		glBindRenderbuffer(GL_RENDERBUFFER, rboColor);
		// Set storage for currently bound renderbuffer.
		glRenderbufferStorage(GL_RENDERBUFFER, GL_RGB, this->window_width, this->window_height);
	}

	// create render buffer for depth
	void createRBOD(void)
	{
		// Depth renderbuffer
		glGenRenderbuffers(1, &rboDepth);
		glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, window_width, window_height);
	}

	// create frame buffer
	void createFBO(void)
	{
		// Framebuffer
		glGenFramebuffers(1, &fbo);
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
		glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rboColor); // Set renderbuffers of color for currently bound framebuffer
		glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth); // Set renderbuffers of depth for currently bound framebuffer

		// Set to write to the framebuffer.
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);

		// Tell glReadPixels where to read from.
		glReadBuffer(GL_COLOR_ATTACHMENT0);
	}

	//Visualize: displays shape points as pixels on GL window
	void visualize(void)
	{
		//------------------------------------ Attach the Shader by Use program ------------------------------------------------
		if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
		{
			//this->setLinghtMaterialR();
			glUseProgram(this->NRGA_shadersR->getProgShaderId());
			glUniform1f(this->NRGA_shadersR->getPositionUniform(), NRGA_sim->NRGA_getParameters()->transparency);
			glUniformMatrix4fv(this->NRGA_shadersR->getProjectionMatrixUniform(), 1, GL_FALSE, this->projection_matrix_);
			glUniformMatrix4fv(this->NRGA_shadersR->getModelViewMatrixUniform(), 1, GL_FALSE, this->modelview_matrix_);
			glUniform4fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Light.Position"), 1, NRGA_sim->NRGA_getParameters()->LightPositionR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Light.La"), 1, NRGA_sim->NRGA_getParameters()->LAmbientColorR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Light.Ld"), 1, NRGA_sim->NRGA_getParameters()->LDeffuseColorR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Light.Ls"), 1, NRGA_sim->NRGA_getParameters()->LSpecularColorR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Material.Ka"), 1, NRGA_sim->NRGA_getParameters()->MAmbientR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Material.Kd"), 1, NRGA_sim->NRGA_getParameters()->MDiffuseR);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Material.Ks"), 1, NRGA_sim->NRGA_getParameters()->MSpecularR);
			glUniform1f(glGetUniformLocation(this->NRGA_shadersR->getProgShaderId(), "Material.Shininess"), NRGA_sim->NRGA_getParameters()->MShinynessR);
		}
		//--------------------------------------- vertex buffer object for reference points ------------------------------------
		glBindVertexArray(this->p_vaoR);
		if (this->NRGA_sim->NRGA_getParameters()->useShadersR == false)
		{
			glBindBuffer(GL_ARRAY_BUFFER, this->p_vboR);
			
			if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
			{
			  glEnableClientState(GL_VERTEX_ARRAY);
			  glVertexPointer(3, GL_DOUBLE, sizeof(this->Rpts3D[0]), NULL);
			  glEnableClientState(GL_COLOR_ARRAY);
			  glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(this->Rpts3D[0]), this->colOffseR);
			}
			else if(this->NRGA_sim->NRGA_getParameters()->dimR == 4)
			{
			  glEnableClientState(GL_VERTEX_ARRAY);
			  glVertexPointer(3, GL_DOUBLE, sizeof(this->reference_points4D[0]), NULL);
			  glEnableClientState(GL_COLOR_ARRAY);
			  glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(this->reference_points4D[0]), this->colOffseR);
			}
			else
			{
			  ;
			}
		}
		if (this->normalsR != NULL && p_numPoints_R > 0)
		{
			glBindBuffer(GL_ARRAY_BUFFER, this->p_nboR);
		}
		if (this->IndicesR != NULL && this->p_numIndices_R > 0)
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboR);
		}
		if (draw_Reference)
		{
			if (toggleCM)
			{
				if (this->p_numIndices_R > 0)
					glDrawElements(GL_TRIANGLES, this->p_numIndices_R, GL_UNSIGNED_INT, 0);
				else
				{
					glDrawArrays(GL_POINTS, 0, p_numPoints_R);
				}

			}
			else
				glDrawArrays(GL_POINTS, 0, p_numPoints_R);
		}
		if (this->NRGA_sim->NRGA_getParameters()->useShadersR == false)
		{
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
		}
		//----------------------------------------Detach the shader in the End of Visualization Loop-----------------------------
		if (this->NRGA_sim->NRGA_getParameters()->useShadersR == true)
		{
			glUseProgram(0);
		}


		//------------------------------------ Attach the Shader by Use program ------------------------------------------------
		if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
		{
			//this->setLinghtMaterialT();
			glUseProgram(this->NRGA_shadersT->getProgShaderId());
			const double timeScale = 0.008f;
			double pT[2] = { 0.5f * sinf(50.0 * timeScale), 0.5f * cosf(50.0 * timeScale) };
			glUniform2dv(this->NRGA_shadersT->getPositionUniform(), 1, (const GLdouble *)&pT);
			glUniformMatrix4fv(this->NRGA_shadersT->getProjectionMatrixUniform(), 1, GL_FALSE, this->projection_matrix_);
			glUniformMatrix4fv(this->NRGA_shadersT->getModelViewMatrixUniform(), 1, GL_FALSE, this->modelview_matrix_);
			glUniform4fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Light.Position"), 1, NRGA_sim->NRGA_getParameters()->LightPositionT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Light.La"), 1, NRGA_sim->NRGA_getParameters()->LAmbientColorT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Light.Ld"), 1, NRGA_sim->NRGA_getParameters()->LDeffuseColorT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Light.Ls"), 1, NRGA_sim->NRGA_getParameters()->LSpecularColorT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Material.Ka"), 1, NRGA_sim->NRGA_getParameters()->MAmbientT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Material.Kd"), 1, NRGA_sim->NRGA_getParameters()->MDiffuseT);
			glUniform3fv(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Material.Ks"), 1, NRGA_sim->NRGA_getParameters()->MSpecularT);
			glUniform1f(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "Material.Shininess"), NRGA_sim->NRGA_getParameters()->MShinynessT);
			glUniform1i(glGetUniformLocation(this->NRGA_shadersT->getProgShaderId(), "ColorMap"), 0);

			glActiveTexture(GL_TEXTURE0 + 0);
			glBindTexture(GL_TEXTURE_1D, this->p_textureId);
		}

		// ------------------------------------vertex buffer object for template points ----------------------------------------
		glBindVertexArray(this->p_vaoT);

		if (this->normalsT != NULL && p_numPoints_T > 0)
		{
			glBindBuffer(GL_ARRAY_BUFFER, this->p_nboT);
		}
		if (this->IndicesT != NULL && this->p_numIndices_T > 0)
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->p_iboT);
		}

		if (draw_Template)
		{
			if (toggleCM)
			{
				if (this->p_numIndices_T > 0)
					glDrawElements(GL_TRIANGLES, this->p_numIndices_T, GL_UNSIGNED_INT, 0);
				else
					glDrawArrays(GL_POINTS, 0, p_numPoints_T);
			}
			else
				glDrawArrays(GL_POINTS, 0, p_numPoints_T);
		}

		//----------------------------------------Detach the shader in the End of Visualization Loop-----------------------------
		if (this->NRGA_sim->NRGA_getParameters()->useShadersT == true)
		{
			glUseProgram(0);
		}
	}

	void drawline(float x1, float y1, float x2, float y2, float z1, float z2)
	{
		glBegin(GL_LINES);
		glVertex3f(x1, y1, z1);
		glVertex3f(x2, y2, z2);
		glEnd();
	}

	void drawgrid()
	{
		glClearColor(0.0, 0.0, 0.0, 1.0);
		glColor3ub(192, 192, 192);

		for (float i = -2; i < 2; i += 0.5f)
		{
		  for(float j = -2; j < 2; j += 0.5f)
		  {
		      glLineWidth(0.1);
		      drawline(-2.0, i, 2.0, i, j, j);
		  }
		}

		for (float i = -2; i < 2; i += 0.5f)
		{
		   for(float j = -2; j < 2; j += 0.5f)
		   {
		    	glLineWidth(0.1);
			drawline(i, -2, i, 2, j, j); 
		   }
		}
	}

	void drawKNN()
    {
        for(int i=0; i < iter+1; i++ ) //p_numPoints_T
        {
            std::vector<size_t> indexesNNT(KDTree3D_T->k);
            std::vector<double> distNNT(KDTree3D_T->k);
            KDTree3D_T->getKNNS(this->Tpts3D[i], indexesNNT, distNNT);

            std::vector<size_t> indexesNNR(KDTree3D_R->k);
            std::vector<double> distNNR(KDTree3D_R->k);     
            KDTree3D_R->getKNNS(this->Tpts3D[i], indexesNNR, distNNR);
            
            double Txmin = this->Tpts3D[i].pos.x;
            double Tymin = this->Tpts3D[i].pos.y; 
            double Tzmin = this->Tpts3D[i].pos.z;
            
            double Txmax = this->Tpts3D[i].pos.x;
            double Tymax = this->Tpts3D[i].pos.y;
            double Tzmax = this->Tpts3D[i].pos.z;
            
            double Rxmin = this->Tpts3D[i].pos.x;
            double Rymin = this->Tpts3D[i].pos.y; 
            double Rzmin = this->Tpts3D[i].pos.z;
            
            double Rxmax = this->Tpts3D[i].pos.x;
            double Rymax = this->Tpts3D[i].pos.y;
            double Rzmax = this->Tpts3D[i].pos.z;
            
            for (int _t = 0; _t < KDTree3D_T->k; _t++)
            {
                if(this->Tpts3D[indexesNNT.at(_t)].pos.x > Txmax)
                    Txmax = this->Tpts3D[indexesNNT.at(_t)].pos.x;
                else if (this->Tpts3D[indexesNNT.at(_t)].pos.x < Txmin)
                    Txmin = this->Tpts3D[indexesNNT.at(_t)].pos.x;
                
                if(this->Tpts3D[indexesNNT.at(_t)].pos.y > Tymax)
                    Tymax = this->Tpts3D[indexesNNT.at(_t)].pos.y;
                else if (this->Tpts3D[indexesNNT.at(_t)].pos.y < Tymin)
                    Tymin = this->Tpts3D[indexesNNT.at(_t)].pos.y;
                
                if(this->Tpts3D[indexesNNT.at(_t)].pos.z > Tzmax)
                    Tzmax = this->Tpts3D[indexesNNT.at(_t)].pos.z;
                else if (this->Tpts3D[indexesNNT.at(_t)].pos.z < Tzmin)
                    Tzmin = this->Tpts3D[indexesNNT.at(_t)].pos.z;
            }
            
            for (int r = 0; r < KDTree3D_R->k; r++)
            {
                if(this->Rpts3D[indexesNNR.at(r)].pos.x > Rxmax)
                    Rxmax = this->Rpts3D[indexesNNR.at(r)].pos.x;
                else if (this->Rpts3D[indexesNNR.at(r)].pos.x < Rxmin)
                    Rxmin = this->Rpts3D[indexesNNR.at(r)].pos.x;
                
                if(this->Rpts3D[indexesNNR.at(r)].pos.y > Rymax)
                    Rymax = this->Rpts3D[indexesNNR.at(r)].pos.y;
                else if (this->Rpts3D[indexesNNR.at(r)].pos.y < Rymin)
                    Rymin = this->Rpts3D[indexesNNR.at(r)].pos.y;
                
                if(this->Rpts3D[indexesNNR.at(r)].pos.z > Rzmax)
                    Rzmax = this->Rpts3D[indexesNNR.at(r)].pos.z;
                else if (this->Rpts3D[indexesNNR.at(r)].pos.z < Rzmin)
                    Rzmin = this->Rpts3D[indexesNNR.at(r)].pos.z;
            }
            
            glColor3ub(12, 218, 225);
            glLineWidth(1.9);
            
            drawline(Txmin, Tymin, Txmax, Tymin, -1.0, -1.0);
            drawline(Txmin, Tymax, Txmax, Tymax, -1.0, -1.0);
            //drawline(Txmin, Tymin, Tzmin,  Txmin, Tymax, Tzmax);
            
            drawline(Txmin, Tymin, Txmin, Tymax, -1.0, -1.0);
            drawline(Txmax, Tymin, Txmax, Tymax, -1.0, -1.0);
            //drawline(Txmax, Tymin, Tzmin,  Txmax, Tymax, Tzmax);
            
            
            
            glColor3ub(255, 0, 0);
            glLineWidth(1.9);
            
            drawline(Rxmin, Rymin, Rxmax, Rymin, -1.0, -1.0);
            drawline(Rxmin, Rymax, Rxmax, Rymax, -1.0, -1.0);
            //drawline(Txmin, Tymin, Tzmin,  Txmin, Tymax, Tzmax);
            
            drawline(Rxmin, Rymin, Rxmin, Rymax, -1.0, -1.0);
            drawline(Rxmax, Rymin, Rxmax, Rymax, -1.0, -1.0);
            //drawline(Txmax, Tymin, Tzmin,  Txmax, Tymax, Tzmax);
            
            //min.x, min.y, min.z
            //min.x, max.y, min.z
            //min.x, min.y, max.z
            //min.x, max.y, max.z
            
            //max.x, min.y, min.z
            //max.x, max.y, min.z
            //max.x, min.y, max.z
            //max.x, max.y, max.z
            
        }
    }
        
    void drawForceField()
    {
        if(this->draw_gravityField == true)
        {                
            this->p_left_bound = -1;//-6;//-1;
            this->p_right_bound = 1;//6;//1;
            this->p_top_bound = 1;//6;//2;
            this->p_bottom_bound = -1;//-6;//-2;
            this->p_nDepth_bound = 1;//2;//2;
            this->p_fDepth_bound = -1;//-2;//1;
            
            double fx           = 0.0;
            double fy           = 0.0;
            double fz           = 0.0;
            double ax           = 0.0;
            double ay           = 0.0;
            double az           = 0.0;
            double vx           = 0.0;
            double vy           = 0.0;
            double vz           = 0.0;
            double posX         = 0.0;
            double posY         = 0.0;
            double posZ         = 0.0;
            double distX        = 0.0;
            double distY        = 0.0;
            double distZ        = 0.0;
            double angle        = 0.0;
            double arrow_length = 0.0000008;//0.08;

            double spacing_X = 0.025;//0.5;
            double spacing_Y = 0.025;//0.5;
            double spacing_Z = 1.0;//0.5;
            
            int nx = (this->p_right_bound - this->p_left_bound)/ spacing_X;
            int ny = (this->p_top_bound - this->p_bottom_bound)/ spacing_Y;
            int nz = (this->p_nDepth_bound - this->p_fDepth_bound) / spacing_Z;

            //glClearColor(0.0, 0.0, 0.0, 1.0);
    glColor3ub(12, 218, 225);
            glLineWidth(0.05);

            glBegin(GL_LINES);
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nz; k++)
                    {
                        fx = 0.0;
                        fy = 0.0;
                        fz = 0.0;
                        
                        for (long r = 0; r < this->p_numPoints_R; r++)
                        {
                            distX = this->p_left_bound +   (i*spacing_X) - this->Rpts3D[r].pos.x; // can be optimized
                            distY = this->p_bottom_bound + (j*spacing_Y) - this->Rpts3D[r].pos.y; // can be optimized
                            distZ = this->p_fDepth_bound + (k*spacing_Z) - this->Rpts3D[r].pos.z; // can be optimized
                            
                            fx += (0.04 * distX) / pow(pow(distX, 2) + pow(distY, 2) + pow(distZ, 2.0) + pow(0.1, 2.0), 1.5);
                            fy += (0.04 * distY) / pow(pow(distX, 2) + pow(distY, 2) + pow(distZ, 2.0) + pow(0.1, 2.0), 1.5);
                            fz += (0.04 * distZ) / pow(pow(distX, 2) + pow(distY, 2) + pow(distZ, 2.0) + pow(0.1, 2.0), 1.5);
                        }

                        vx = -fx * 0.0008;
                        vy = -fy * 0.0008;
                        vz = -fz * 0.0008;
                        
                        double a = std::sqrt(pow(fx,2.0) + pow(fy,2.0) + pow(fz,2.0));
                        posX = (this->p_left_bound   + (i*spacing_X)) + (vx * 3.0);
                        posY = (this->p_bottom_bound + (j*spacing_Y)) + (vy * 3.0);
                        posZ = (this->p_fDepth_bound + (k*spacing_Z)) + (vz * 3.0);
                        
                        distX = posX - (this->p_left_bound + (i*spacing_X));
                        distY = posY - (this->p_bottom_bound + (j*spacing_Y));
                        distZ = posZ - (this->p_fDepth_bound + (k*spacing_Z));

                        angle = atan2(distY, distX);

                        glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
                        glVertex3d((GLdouble)(this->p_left_bound + (i*spacing_X)), (GLdouble)(this->p_bottom_bound + (j*spacing_Y)), (GLdouble)(this->p_fDepth_bound + (k*spacing_Z)));

//                             // angle for 1st quadrant
//                             if ((posX > 0.0 && posY > 0.0) && (distX <= 0.0 && distY <= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle + 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle - 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY > 0.0) && (distX > 0.0 && distY > 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY > 0.0) && (distX > 0.0 && distY < 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY > 0.0) && (distX < 0.0 && distY > 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
// 
//                             // angle for 2nd quadrant
//                             if ((posX < 0.0 && posY > 0.0) && (distX >= 0.0 && distY <= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY > 0.0) && (distX >= 0.0 && distY >= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY > 0.0) && (distX <= 0.0 && distY <= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle + 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle - 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY > 0.0) && (distX <= 0.0 && distY >= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
// 
//                             // angle for 3rd quadrant
//                             if ((posX < 0.0 && posY < 0.0) && (distX > 0.0 && distY > 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY < 0.0) && (distX <= 0.0 && distY <= 0.0))
//                             {
//                                
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle + 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle - 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY < 0.0) && (distX >= 0.0 && distY <= 0.0))
//                             {
//                                
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX < 0.0 && posY < 0.0) && (distX <= 0.0 && distY >= 0.0))
//                             {
//                                
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
// 
//                             // angle for 4th quadrant
//                             if ((posX > 0.0 && posY < 0.0) && (distX < 0.0 && distY > 0.0))
//                             {
//                             
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY < 0.0) && (distX >= 0.0 && distY <= 0.0))
//                             {
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)),posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY < 0.0) && (distX >= 0.0 && distY >= 0.0))
//                             {
//                                 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle + 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX + arrow_length * cos(3.14 + angle - 3.14 / 9)), (GLdouble)(posY + arrow_length * sin(3.14 + angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
//                             else if ((posX > 0.0 && posY < 0.0) && (distX <= 0.0 && distY <= 0.0))
//                             {
//                                 
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle + 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle + 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
// 
//                                 glVertex3d((GLdouble)(posX - arrow_length * cos(angle - 3.14 / 9)), (GLdouble)(posY - arrow_length * sin(angle - 3.14 / 9)), posZ);
//                                 glVertex3d((GLdouble)posX, (GLdouble)posY, (GLdouble)posZ);
//                             }
                    }
                }
            }
            glEnd();
        }
    }
	
	bool map_to_sphere(const QPoint& _v2D, Eigen::Vector3d& _v3D)
	{
		// This is actually doing the Sphere/Hyperbolic sheet hybrid thing,
		// based on Ken Shoemake's ArcBall in Graphics Gems IV, 1993.
		double x = (2.0*_v2D.x() - width()) / width();
		double y = -(2.0*_v2D.y() - height()) / height();
		double xval = x;
		double yval = y;
		double x2y2 = xval*xval + yval*yval;

		const double rsqr = this->TRACKBALL_RADIUS * this->TRACKBALL_RADIUS;
		_v3D[0] = xval;
		_v3D[1] = yval;
		if (x2y2 < 0.5*rsqr) {
			_v3D[2] = sqrt(rsqr - x2y2);
		}
		else {
			_v3D[2] = 0.5*rsqr / sqrt(x2y2);
		}

		return true;
	}

	void set_scene_pos(const Eigen::Vector3d& _cog, float _radius)
	{
		center_ = _cog;
		radius_ = _radius;
		glFogf(GL_FOG_START, 1.5*_radius);
		glFogf(GL_FOG_END, 3.0*_radius);

		update_projection_matrix();
		view_all();
	}

	void update_projection_matrix()
	{
		makeCurrent();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(this->FOV, (GLfloat)this->window_width / (GLfloat)this->window_height, -100.1f * this->radius_, (this->window_width / 10.0f) * this->radius_);
		glGetFloatv(GL_PROJECTION_MATRIX, this->projection_matrix_);
		glMatrixMode(GL_MODELVIEW);
	}

	void translate(const Eigen::Vector3d& _trans)
	{
		// Translate the object by _trans
		// Update modelview_matrix_
		makeCurrent();
		glLoadIdentity();
		glTranslated(_trans[0], _trans[1], _trans[2]);
		glMultMatrixf(this->modelview_matrix_);
		glGetFloatv(GL_MODELVIEW_MATRIX, this->modelview_matrix_);
	}

	void rotate(const Eigen::Vector3d& _axis, float _angle)
	{
		// Rotate around center center_, axis _axis, by angle _angle
		// Update modelview_matrix_

		Eigen::Vector3d t(
			modelview_matrix_[0] * this->center_[0] +
			modelview_matrix_[4] * this->center_[1] +
			modelview_matrix_[8] * this->center_[2] +
			modelview_matrix_[12],

			modelview_matrix_[1] * this->center_[0] +
			modelview_matrix_[5] * this->center_[1] +
			modelview_matrix_[9] * this->center_[2] +
			modelview_matrix_[13],

			modelview_matrix_[2] * this->center_[0] +
			modelview_matrix_[6] * this->center_[1] +
			modelview_matrix_[10] * this->center_[2] +
			modelview_matrix_[14]);

		makeCurrent();
		glLoadIdentity();
		
		glTranslatef(t[0], t[1], t[2]);
		glRotated(_angle, _axis[0], _axis[1], _axis[2]);
		glTranslatef(-t[0], -t[1], -t[2]);
		glMultMatrixf(modelview_matrix_);
		glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix_);
	}

	void view_all()
	{
		translate(Eigen::Vector3d(
			-(modelview_matrix_[0] * center_[0] +
			modelview_matrix_[4] * center_[1] +
			modelview_matrix_[8] * center_[2] +
			modelview_matrix_[12]),

			-(modelview_matrix_[1] * center_[0] +
			modelview_matrix_[5] * center_[1] +
			modelview_matrix_[9] * center_[2] +
			modelview_matrix_[13]),

			-(modelview_matrix_[2] * center_[0] +
			modelview_matrix_[6] * center_[1] +
			modelview_matrix_[10] * center_[2] +
			modelview_matrix_[14] +
			3.0*radius_)));
	}

	// save image Sequence 
	void slotSnapshot(void)
	{
		QImage image;
		size_t w(width()), h(height());
		GLenum buffer(GL_BACK);

		try
		{
			image = QImage(w, h, QImage::Format_RGB32);
			std::vector<GLubyte> fbuffer(3 * w*h);

			qApp->processEvents();
			makeCurrent();
			updateGL();
			glFinish();

			glReadBuffer(buffer);
			glPixelStorei(GL_PACK_ALIGNMENT, 1);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			paintGL();
			glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &fbuffer[0]);

			unsigned int x, y, offset;
			for (y = 0; y < h; ++y) {
				for (x = 0; x < w; ++x) {
					offset = 3 * (y*w + x);
					image.setPixel(x, h - y - 1, qRgb(fbuffer[offset],
						fbuffer[offset + 1],
						fbuffer[offset + 2]));
				}
			}


			QString name = QString(this->NRGA_sim->NRGA_getParameters()->out_registered.c_str()) + "snapshot-";
#if defined(_MSC_VER)
			{
				std::stringstream s;
				QDateTime         dt = QDateTime::currentDateTime();
				s << dt.date().year()
					<< std::setw(2) << std::setfill('0') << dt.date().month()
					<< std::setw(2) << std::setfill('0') << dt.date().day()
					<< std::setw(2) << std::setfill('0') << dt.time().hour()
					<< std::setw(2) << std::setfill('0') << dt.time().minute()
					<< std::setw(2) << std::setfill('0') << dt.time().second();
				name += QString(s.str().c_str());
			}
#else
			//name += QDateTime::currentDateTime().toString("yyMMddhhmmss");
			if(iter < 10)
			    name += QString("000") + QString::number(iter);
			else if(iter >= 10 && iter < 100)
			    name += QString("00") + QString::number(iter);
			else if(iter >= 100 && iter < 1000)
			    name += QString("0") + QString::number(iter);
#endif
			name += ".png";

			image.save(name, "PNG");
			}
		catch (std::bad_alloc&)
		{
			qWarning("Mem Alloc Error");
		}

		}

	// save template ply file
	void saveTemplate(void)
	{
		try
		{
			qApp->processEvents();
			makeCurrent();

			if(this->NRGA_sim->NRGA_getParameters()->dimT == 3)
			{
			    Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				&this->Tpts3D[0].pos.x,
				p_numPoints_T,
				3,
				Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
				);
			    NRGA_plyModuleT->setVertices(T_pts);
                            
                            if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 0)
                            {
                                Eigen::MatrixXi Same_color(p_numPoints_T, 3);
                                Eigen::VectorXi color_(3); color_ << (iter), 0, (255 - (iter));
                                Same_color.setZero();
                                Same_color.rowwise() += color_.transpose();
                                NRGA_plyModuleT->setColors(Same_color);
                            }
                            else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 1)
                            {
                                Eigen::MatrixXi Same_color = (Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                                    &this->Tpts3D[0].color.R,
                                    p_numPoints_T,
                                    3,
                                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                                    )).cast<int>();
                                NRGA_plyModuleT->setColors(Same_color);
                            }
                            else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 2)
                            {
                                ;
                            } 
			}
			else if(this->NRGA_sim->NRGA_getParameters()->dimT == 4)
			{
			    Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				&this->template_points4D[0].pos.x,
				p_numPoints_T,
				3,
				Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
				);
			    NRGA_plyModuleT->setVertices(T_pts);
                            
                            if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 0)
                            {
                                Eigen::MatrixXi Same_color(p_numPoints_T, 3);
                                Eigen::VectorXi color_(3); color_ << (iter), 0, (255 - (iter));
                                Same_color.setZero();
                                Same_color.rowwise() += color_.transpose();
                                NRGA_plyModuleT->setColors(Same_color);
                            }
                            else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 1)
                            {
                                Eigen::MatrixXi Same_color = (Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                                    &this->template_points4D[0].color.R,
                                    p_numPoints_T,
                                    3,
                                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
                                    )).cast<int>();
                                NRGA_plyModuleT->setColors(Same_color);
                            }
                            else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 2)
                            {
                                ;
                            } 
			}
                        
			QString name = QString(this->NRGA_sim->NRGA_getParameters()->out_registered.c_str()) 
				     + "Template-" + QString::number(iter) + "_";

#if defined(_MSC_VER)
			{
			    std::stringstream s;
			    QDateTime         dt = QDateTime::currentDateTime();
			    s << dt.date().year()
				    << std::setw(2) << std::setfill('0') << dt.date().month()
				    << std::setw(2) << std::setfill('0') << dt.date().day()
				    << std::setw(2) << std::setfill('0') << dt.time().hour()
				    << std::setw(2) << std::setfill('0') << dt.time().minute()
				    << std::setw(2) << std::setfill('0') << dt.time().second();
			    name += QString(s.str().c_str());
			}
#else
			//name += QDateTime::currentDateTime().toString("yyMMddhhmmss");
			if(iter < 10)
			    name += QString("000") + QString::number(iter);
			else if(iter >= 10 && iter < 100)
			    name += QString("00") + QString::number(iter);
			else if(iter >= 100 && iter < 1000)
			    name += QString("0") + QString::number(iter);
#endif
			name += ".ply";
			bool successT = NRGA_plyModuleT->writePLY(name.toStdString(), true, true, false, false, false);
			
			if (successT == true)
				std::cout << "Template file has been read properly" << std::endl;
			
                        // For The first iteration we save the "Reference" too
			if(iter == 1)
			{	  
			    if(this->NRGA_sim->NRGA_getParameters()->dimR == 3)
			    {
				Eigen::MatrixXd R_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				    &this->Rpts3D[0].pos.x,
				    p_numPoints_R,
				    3,
				    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
				    );
				NRGA_plyModuleR->setVertices(R_pts);
                                
                                if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 0)
                                {
                                    Eigen::MatrixXi Same_color(p_numPoints_R, 3);
                                    Eigen::VectorXi color_(3); color_ << (iter), 0, (255 - (iter));
                                    Same_color.setZero();
                                    Same_color.rowwise() += color_.transpose();
                                    NRGA_plyModuleR->setColors(Same_color);
                                }
                                else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 1)
                                {
                                    Eigen::MatrixXi Same_color = (Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                                        &this->Rpts3D[0].color.R,
                                        p_numPoints_R,
                                        3,
                                        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                                        )).cast<int>();
                                    NRGA_plyModuleR->setColors(Same_color);
                                }
                                else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 2)
                                {
                                    ;
                                } 
			    }
			    else if(this->NRGA_sim->NRGA_getParameters()->dimR == 4)
			    {
				Eigen::MatrixXd R_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
				    &this->reference_points4D[0].pos.x,
				    p_numPoints_R,
				    3,
				    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
				    );
				NRGA_plyModuleR->setVertices(R_pts);
                                if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 0)
                                {
                                    Eigen::MatrixXi Same_color(p_numPoints_R, 3);
                                    Eigen::VectorXi color_(3); color_ << (iter), 0, (255 - (iter));
                                    Same_color.setZero();
                                    Same_color.rowwise() += color_.transpose();
                                    NRGA_plyModuleR->setColors(Same_color);
                                }
                                else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 1)
                                {
                                    Eigen::MatrixXi Same_color = (Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                                        &this->reference_points4D[0].color.R,
                                        p_numPoints_R,
                                        3,
                                        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                                        )).cast<int>();
                                    NRGA_plyModuleR->setColors(Same_color);
                                }
                                else if(this->NRGA_sim->NRGA_getParameters()->p_color_output == 2)
                                {
                                    ;
                                } 
			    }

			    QString name = QString(this->NRGA_sim->NRGA_getParameters()->out_registered.c_str()) 
				      + "Reference-" + QString::number(iter);
#if defined(_MSC_VER)
			    {
				std::stringstream s;
				QDateTime         dt = QDateTime::currentDateTime();
				s << dt.date().year()
					<< std::setw(2) << std::setfill('0') << dt.date().month()
					<< std::setw(2) << std::setfill('0') << dt.date().day()
					<< std::setw(2) << std::setfill('0') << dt.time().hour()
					<< std::setw(2) << std::setfill('0') << dt.time().minute()
					<< std::setw(2) << std::setfill('0') << dt.time().second();
				name += QString(s.str().c_str());
			    }
#else
			    name += QDateTime::currentDateTime().toString("yyMMddhhmmss");
#endif
			    name += ".ply";
			    bool successR = NRGA_plyModuleR->writePLY(name.toStdString(), true, true, false, false, false);
			    if (successR == true)
			      std::cout << "Reference file has been read properly" << std::endl;
			}
			
		}
		catch (std::bad_alloc&)
		{
			qWarning("Mem Alloc Error");
		}

	}

	// key press event
	void keyPressEvent(QKeyEvent *keyEvent)
	{       
            if (keyEvent->type() == QEvent::KeyPress)
            {
                switch (keyEvent->key())
                {
                    case Qt::Key_R:
                            this->R_Color = !this->R_Color;
                            this->updateGL();
                            break;
                    case Qt::Key_G:
                            this->G_Color = !this->G_Color;
                            updateGL();
                            break;
                    case Qt::Key_B:
                            this->B_Color = !this->B_Color;
                            updateGL();
                            break;
                    case Qt::Key_A:
                            this->draw_axis = !this->draw_axis;
                            std::cout << "A pressed to draw Axes" << std::endl;
                            updateGL();
                            break;
                    case Qt::Key_V:
                            this->draw_grid = !this->draw_grid;
                            std::cout << "key pressed" << std::endl;
                            updateGL();
                            break;
                     case Qt::Key_K:
                            this->draw_KNN = !this->draw_KNN;
                            std::cout << "K pressed to show KNN" << std::endl;
                            updateGL();
                            break;
                    case Qt::Key_I:
                            this->draw_Image = !this->draw_Image;
                            std::cout << "I pressed to save images" << std::endl;
                            break;
                    case Qt::Key_O:
                            this->draw_shape = !this->draw_shape;
                            std::cout << "O pressed to draw objects" << std::endl;
                            updateGL();
                            break;
                    case Qt::Key_P:
                            this->save_ply = !this->save_ply;
                            std::cout << "P pressed to save ply" << std::endl;
                            updateGL();
                            break;
                    case Qt::Key_F:
                            this->draw_gravityField = !this->draw_gravityField;
                            std::cout << "G pressed to draw gravitational vectors" << std::endl;
                            updateGL();
                            break;
                    case Qt::Key_T:
                            this->draw_Reference = !this->draw_Reference;
                            updateGL();
                            break;
                    case Qt::Key_Space:
                            this->draw_Simulation = !this->draw_Simulation;
                            updateGL();
                            break;
                    case Qt::Key_S:
                            this->draw_Template = !this->draw_Template;
                            updateGL();
                            break;
                    case Qt::Key_C:
                            this->toggleCM = !this->toggleCM;
                            updateGL();
                            break;
                    case Qt::Key_H:
                            std::cout << "Keys:\n";
                            std::cout << "  I\t Make snapshot\n";
                            std::cout << "  C\t enable/disable back face culling\n";
                            std::cout << "  F\t enable/disable fog\n";
                            std::cout << "  H\t Display information\n";
                            std::cout << "  N\t enable/disable display of vertex normals\n";
                            std::cout << "  T\t Toggle cloud and mesh visualization\n";
                            break;
                    case Qt::Key_Q:
                    case Qt::Key_Escape:
                            qApp->quit();
                }

                keyEvent->ignore();
            }
	}

	void keyReleaseEvent(QKeyEvent *keyEvent)
    {
        if (keyEvent->type() == QEvent::KeyPress)
        {
            switch (keyEvent->key())
            {
                case Qt::Key_R:
                        this->R_Color = !this->R_Color;
                        updateGL();
                        break;
                case Qt::Key_G:
                        this->G_Color = this->G_Color;
                        updateGL();
                        break;
                case Qt::Key_B:
                        this->B_Color = this->B_Color;
                        updateGL();
                        break;
            }
            
            keyEvent->ignore();
        }
    }

	// mouse press event
	void mousePressEvent(QMouseEvent *mouseEvent)
	{
		// popup menu
		if (mouseEvent->button() == Qt::RightButton && mouseEvent->buttons() == Qt::RightButton)
		{
			//popup_menu_->exec(QCursor::pos());
		}
		else
		{
			last_point_ok_ = map_to_sphere(last_point_2D_ = mouseEvent->pos(), last_point_3D_);
		}
		mouseEvent->accept();
		updateGL();
	}

	// mouse release event
	void mouseReleaseEvent(QMouseEvent *mouseEvent)
	{
		this->last_point_ok_ = false;

		mouseEvent->accept();
		updateGL();
	}

	// mouse move event
	void mouseMoveEvent(QMouseEvent *mouseEvent)
	{
		QPoint newPoint2D = mouseEvent->pos();

		// Left button: rotate around center_
		// Middle button: translate object
		// Left & middle button: zoom in/out


		Eigen::Vector3d  newPoint3D;
		bool   newPoint_hitSphere = map_to_sphere(newPoint2D, newPoint3D);

		float dx = newPoint2D.x() - this->last_point_2D_.x();
		float dy = newPoint2D.y() - this->last_point_2D_.y();

		float w = width();
		float h = height();

		// enable GL context
		makeCurrent();

		// move in z direction
		if ((mouseEvent->buttons() == (Qt::LeftButton + Qt::MidButton)) ||
			(mouseEvent->buttons() == Qt::LeftButton && mouseEvent->modifiers() == Qt::ControlModifier))
		{
			float value_y = radius_ * dy * w / h;
			translate(Eigen::Vector3d(0.0, 0.0, value_y));
		}


		// move in x,y direction
		else if ((mouseEvent->buttons() == Qt::MidButton) ||
			(mouseEvent->buttons() == Qt::LeftButton && mouseEvent->modifiers() == Qt::AltModifier))
		{
			float z = -(
				modelview_matrix_[2] * center_[0] +
				modelview_matrix_[6] * center_[1] +
				modelview_matrix_[10] * center_[2] +
				modelview_matrix_[14]) /
				(modelview_matrix_[3] * center_[0] +
				modelview_matrix_[7] * center_[1] +
				modelview_matrix_[11] * center_[2] +
				modelview_matrix_[15]
				);

			float aspect = w / h;
			float near_plane = 0.01 * radius_;
			float top = tan(this->FOV / 2.0f*M_PI / 180.0f) * near_plane;
			float right = aspect*top;

			translate(Eigen::Vector3d(2.0*dx / w*right / near_plane*z,
				-2.0*dy / h*top / near_plane*z,
				0.0f));
		}

		// rotate
		else if (mouseEvent->buttons() == Qt::LeftButton) {

			if (last_point_ok_) {
				if ((newPoint_hitSphere = map_to_sphere(newPoint2D, newPoint3D))) {

					Eigen::Vector3d axis(last_point_3D_ - newPoint3D);
					if (axis.norm() < 1e-4) {
						axis = Eigen::Vector3d(1, 0, 0);
					}
					else {
						axis.normalize();
					}
					// find the amount of rotation
					Eigen::Vector3d d = last_point_3D_ - newPoint3D;
					float t = 0.5 * d.norm() / TRACKBALL_RADIUS;
					if (t < -1.0)
						t = -1.0;
					else if (t > 1.0)
						t = 1.0;
					float phi = 2.0 * asin(t);
					float angle = phi * 180.0 / M_PI;
					rotate(axis, angle);
				}
			}

		}


		// remember this point4D
		last_point_2D_ = newPoint2D;
		last_point_3D_ = newPoint3D;
		last_point_ok_ = newPoint_hitSphere;

		// trigger redraw
		updateGL();
	}

	//mouse wheel event
	void wheelEvent(QWheelEvent *mouseEvent)
	{
            if(mouseEvent->modifiers().testFlag(Qt::ControlModifier))
            {
                if (mouseEvent->delta() > 0)  
                {
                    if(R_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[0] += 0.01;
                    if(G_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[1] += 0.01;
                    if(B_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[2] += 0.01;
                    
                     if(R_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[0] += 0.01;
                    if(G_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[1] += 0.01;
                    if(B_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[2] += 0.01;
                }
                else
                {
                    if(R_Color == true  && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[0] -= 0.01;
                    if(G_Color == true  && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[1] -= 0.01;
                    if(B_Color == true  && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientR[2] -= 0.01;
                    
                    if(R_Color == true  && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[0] -= 0.01;
                    if(G_Color == true  && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[1] -= 0.01;
                    if(B_Color == true  && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MAmbientT[2] -= 0.01;
                }
        
            }
            if(mouseEvent->modifiers().testFlag(Qt::ShiftModifier))
            {
                if (mouseEvent->delta() > 0)  
                {
                     if(R_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[0] += 0.01;
                     if(G_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[1] += 0.01;
                     if(B_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[2] += 0.01;
                     
                      if(R_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[0] += 0.01;
                     if(G_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[1] += 0.01;
                     if(B_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[2] += 0.01;
                }
                else
                {
                    if(R_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[0] -= 0.01;
                    if(G_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[1] -= 0.01;
                    if(B_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseR[2] -= 0.01;
                    
                    if(R_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[0] -= 0.01;
                    if(G_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[1] -= 0.01;
                    if(B_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MDiffuseT[2] -= 0.01;
                }
        
            }
            if(mouseEvent->modifiers().testFlag(Qt::AltModifier))
            {
                if (mouseEvent->delta() > 0)  
                {
                     if(R_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[0] += 0.01;
                     if(G_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[1] += 0.01;
                     if(B_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[2] += 0.01;
                     
                      if(R_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[0] += 0.01;
                     if(G_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[1] += 0.01;
                     if(B_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[2] += 0.01;
                }
                else
                {
                    if(R_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[0] -= 0.01;
                    if(G_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[1] -= 0.01;
                    if(B_Color == true && this->draw_Reference == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularR[2] -= 0.01;
                    
                    if(R_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[0] -= 0.01;
                    if(G_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[1] -= 0.01;
                    if(B_Color == true && this->draw_Template == true)
                        this->NRGA_sim->NRGA_getParameters()->MSpecularT[2] -= 0.01;
                }
        
            }
            updateGL();
	}
};


#ifdef _WIN32
int main(int argc, char *argv[]){

	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	SetConsoleTextAttribute(hConsole, 2);
	std::cout << std::endl; std::cout << std::endl;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||*************Non-Rigid Gravitational Approach for Point Set Registration************||" << std::endl;
	std::cout << "||***************Sk Aziz Ali, Vladislav Golyanik, Didier Stricker*********************||" << std::endl;
	std::cout << "||************{Sk_Aziz.Ali, Vladislav.Golyanik, Didier.Stricker}@dfki.de;*************||" << std::endl;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||-----This application is being implemented to solve Non Rigid Registration using----||" << std::endl;
	std::cout << "||---Physically based approach, Gravitational Approach (Non Accelerated version ).----||" << std::endl;
	std::cout << "||-Any modification or usage without permission of the authos violates copy-right law-||" << std::endl;
	std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;

	SetConsoleTextAttribute(hConsole, 5);
	char type;
	std::cout << " For GUI version: ENTER \"G\" , For CONSOLE version: ENTER \"C\" :  ";
	std::cin >> type;

	if (type == 'G' || type == 'g')
	{
		const char* S_InFile1 = "C:\\Users\\saali14\\Desktop\\NRGA\\DATA\\INPUT\\face_t_subsamp.ply";
		const char* S_InFile2 = "C:\\Users\\saali14\\Desktop\\NRGA\\DATA\\INPUT\\face_r_subsamp.ply";
		const char* S_outFile1 = "/home/ali/PHD/NRGA_Application/DATA/OUTPUT/result.ply";

		SetConsoleTextAttribute(hConsole, 100);
		ply_module* _plyModuleT = new ply_module();
		ply_module* _plyModuleR = new ply_module();
		_plyModuleT->readPLY(S_InFile1, true, true, true);
		_plyModuleR->readPLY(S_InFile2, true, true, true);

		GA<double>* ga = new GA<double>();
		NRGA<double>* nrga = new NRGA<double>();

		shaders* shaderR = new shaders();
		shaders* shaderT = new shaders();

		parameters<double>* param = new parameters<double>();
		param->p_G = 6.67 * std::pow(10.0f, -1.0); // -13.0
		param->p_epsilon = 0.2f;
		param->p_dilatation = 0.20f;
		param->p_timeStep = 0.004;  // 0.004
		param->p_spring_Constraint = 0.00008;
		param->p_damping_Constraint = 0.1;
		param->p_bend_Constraint = 2.0f;
		param->k_T = 5;
		param->k_R = 5;
		param->dimT = 3;
		param->dimR = 3;
		param->useShadersT = true;
		param->useShadersR = true;
		param->shaderPath_TMeshVert = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\vertex_shader_Template.glsl";
		param->shaderPath_RMeshVert = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\vertex_shader_Reference.glsl";
		param->shaderPath_TPointsVert = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\vertex_shader_TemplatePoints.glsl";
		param->shaderPath_RPointsVert = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\vertex_shader_ReferencePoints.glsl";
		param->shaderPath_TMeshFrag = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\fragment_shader_Template.glsl";
		param->shaderPath_RMeshFrag = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\fragment_shader_Reference.glsl";
		param->shaderPath_TPointsFrag = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\fragment_shader_TemplatePoints.glsl";
		param->shaderPath_RPointsFrag = "C:\\Users\\saali14\\Desktop\\NRGA\\shaders\\fragment_shader_ReferencePoints.glsl";
		param->bgColor[0] = 0.0f;
		param->bgColor[1] = 0.0f;
		param->bgColor[2] = 0.0f;
		param->bgColor[3] = 0.0f;
		param->pointSize = 5;

		ga->GA_setParameters(param);
		nrga->NRGA_setParameters(param);

		QApplication app(argc, argv);
		GLWidget* widget = new GLWidget(nrga, ga, shaderR, shaderT, _plyModuleT, _plyModuleR);
		widget->resize(640, 900);
		widget->show();

		return app.exec();
	}
	else if (type == 'C' || type == 'c')
	{
		return 0;
	}
	else
	{
		;
	}
}
#endif

void appendStatistics(std::string path, std::string mean , std::string RMSE)
{
   std::ofstream statistics;
   statistics.open(path.c_str(), std::ios_base::app);
   statistics << mean << " " << RMSE << std::endl;
   statistics.close();
}


#ifdef __linux__
int main(int argc, char *argv[]){
	
	std::cout << std::endl; std::cout << std::endl << ANSI_COLOR_YELLOW;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||*************Non-Rigid Gravitational Approach for Point Set Registration************||" << std::endl;
	std::cout << "||***************Sk Aziz Ali, Vladislav Golyanik, Didier Stricker*********************||" << std::endl;
	std::cout <<                                    ANSI_COLOR_GREEN;
	std::cout << "||*******saali14@rhrk.uni-kl.de;{Vladislav.Golyanik, Didier.Stricker}@dfki.de;********||" << std::endl;
	std::cout <<                                    ANSI_COLOR_YELLOW;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||************************************************************************************||" << std::endl;
	std::cout << "||-----This application is being implemented to solve Non Rigid Registration using----||" << std::endl;
	std::cout << "||---Physically based approach, Gravitational Approach (Non Accelerated version ).----||" << std::endl;
	std::cout << "|| Any modification or usage without permission of the authors violates copy-right law||" << std::endl;
	std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
	std::cout <<                                    ANSI_COLOR_RESET;

	char type;
	std::cout << " For GUI version: ENTER \"G\" , For CONSOLE version: ENTER \"C\" :  ";
	std::cin >> type;

	if (type == 'G' || type == 'g')
	{
        parameters<double>* NRGA_param = new parameters<double>();
        NRGA_param->load_settings(std::string(argv[1]));

        parameters<double>* GA_param = new parameters<double>();
        GA_param->load_settings(std::string(argv[2]));
	      
	    GA<double>* 	_GA 	= new GA<double>();
	    NRGA<double>* 	_NRGA 	= new NRGA<double>();

	    shaders* shaderR = new shaders();
	    shaders* shaderT = new shaders();

	    _GA->GA_setParameters(GA_param);
	    _NRGA->NRGA_setParameters(NRGA_param);

	    ply_module* _plyModuleT = new ply_module();
	    ply_module* _plyModuleR = new ply_module();
	    _plyModuleT->readPLY(NRGA_param->in_template.c_str(), true, true, true, false, true);
	    _plyModuleR->readPLY(NRGA_param->in_reference.c_str(), true, true, true, false, true);

	    QApplication app(argc, argv);
	    GLWidget* widget = new GLWidget(_NRGA, _GA, shaderR, shaderT, _plyModuleT, _plyModuleR);
	    widget->resize(640, 1000);
	    widget->show();

	    return app.exec();
	}
	else if (type == 'C' || type == 'c')
	{        
        parameters<double>* NRGA_param = new parameters<double>();
        NRGA_param->load_settings(std::string(argv[1]));

        parameters<double>* GA_param = new parameters<double>();
        GA_param->load_settings(std::string(argv[2]));

        GA<double>* 	_GA 	= new GA<double>();
        NRGA<double>* 	_NRGA 	= new NRGA<double>();
	    
        int reps = 0;
        do
        {
            if(reps++ >= NRGA_param->p_Reps)
                break;
            
            int NUM_PTS_T = 0;
            int NUM_PTS_R = 0;
            int iter = 0;
            kdt_nanoFlann4D<double>* KDTree4D_T = NULL;
            kdt_nanoFlann4D<double>* KDTree4D_R = NULL;
            kdt_nanoFlann<double>*   KDTree3D_T = NULL;
            kdt_nanoFlann<double>*   KDTree3D_R = NULL;
            
            shaders* shaderR = new shaders();
            shaders* shaderT = new shaders();

            _GA->GA_setParameters(GA_param);
            _NRGA->NRGA_setParameters(NRGA_param);

            ply_module* _plyModuleT = new ply_module();
            ply_module* _plyModuleR = new ply_module();
            _plyModuleT->readPLY(NRGA_param->in_template.c_str(), true, true, true, false, true);
            _plyModuleR->readPLY(NRGA_param->in_reference.c_str(), true, true, true, false, true);
            
            
            // shape points color
            p_color colorT = { (unsigned char)0, (unsigned char)255, (unsigned char)0, (unsigned char)1 };
            p_color colorR = { (unsigned char)255, (unsigned char)0, (unsigned char)0, (unsigned char)1 };
            
            velocity<double> vel;
            vel.vx = 0.0f;
            vel.vy = 0.0f;
            vel.vz = 0.0f;

            acceleration<double>  accl;
            accl.ax = 0.0f;
            accl.ay = 0.0f;
            accl.az = 0.0f;

            double p_mass = 1.0f;

            shape_point<double> template_pts;
            shape_point<double> reference_pts;

            std::string path_Reps =   NRGA_param->out_registered + "/" 
                                    + std::to_string(reps) + "_"
                                    + std::to_string((int)NRGA_param->p_noise_percentage)
                                    + NRGA_param->p_noise_types;	   	  
            int createDir2 =   mkdir(path_Reps.c_str(), 0777);

            // Setup the Template Y
            if(NRGA_param->dimT == 3)
            {
                if(_plyModuleT->getColors().size() == 0)
                    template_pts = shape_point<double>(_plyModuleT->getVertices(), vel, accl, colorT, p_mass);
                else
                    template_pts = shape_point<double>(_plyModuleT->getVertices(), vel, accl, _plyModuleT->getColors(), p_mass);
                
                template_pts.applyMissingChunk(25);
                
                NUM_PTS_T  = template_pts.pointSets3D.size();
                
                if(NRGA_param->p_normalize_flag == 1)
                    template_pts.normalizePoints();
                
                
                unsigned seedR = std::chrono::system_clock::now().time_since_epoch().count();
                std::mt19937 generatorR(seedR);
                std::uniform_real_distribution<double> uniformDistRx(-M_PI * (NRGA_param->p_rotation_max / 180.0) , M_PI * (NRGA_param->p_rotation_max / 180.0));
                std::uniform_real_distribution<double> uniformDistRy(-M_PI * (NRGA_param->p_rotation_max / 180.0) , M_PI * (NRGA_param->p_rotation_max / 180.0));
                std::uniform_real_distribution<double> uniformDistRz(-M_PI * (NRGA_param->p_rotation_max / 180.0) , M_PI * (NRGA_param->p_rotation_max / 180.0));
                
                
                Eigen::AngleAxisd rollAngle(uniformDistRx(generatorR),  Eigen::Vector3d::UnitX());
                Eigen::AngleAxisd pitchAngle(uniformDistRy(generatorR), Eigen::Vector3d::UnitY());
                Eigen::AngleAxisd yawAngle(uniformDistRz(generatorR),  Eigen::Vector3d::UnitZ());
                Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle;
                Eigen::Matrix3d rotationMatrix = q.toRotationMatrix();
                template_pts.rotateAbout_XYZ_axis(rotationMatrix);

                if(NRGA_param->p_addnoise_flag == 1)
                {
                    if(NRGA_param->p_noise_types == "G")
                    {
                        template_pts.add3DGaussianNoise(NRGA_param->p_noiseG_variance, 
                                                        template_pts.calcNumberofNoisyPoints(NRGA_param->p_noise_percentage,template_pts.pointSets3D.size()), 
                                                        vel, accl, colorT, p_mass);
                    }
                    if(NRGA_param->p_noise_types == "U")
                    {
                        template_pts.addUniformNoise3D(template_pts.calcNumberofNoisyPoints(NRGA_param->p_noise_percentage,template_pts.pointSets3D.size()), 
                                                        vel, accl, colorT, p_mass);
                    }
                    if(NRGA_param->p_noise_types == "P")
                    {
                        double step_size = 0.000005; 
                        double max_deviation = reps * step_size;
                        template_pts.addPerturbations(max_deviation);
                    }
                }

                _NRGA->NRGA_initShapePoints_T3D(template_pts);
                
            }
            else if(NRGA_param->dimT == 4)
            {
                if(_plyModuleT->getColors().size() > 0 && _plyModuleT->getCurvatures().size() > 0)
                    template_pts = shape_point<double>(_plyModuleT->getVertices(), _plyModuleT->getCurvatures(), vel, accl, _plyModuleT->getColors(), p_mass);
                else
                    std::cout << "Vertex Attributes is NOT initialized for TEMPLATE" << std::endl;
            }
            
            
            // Setup the Reference X
            if(NRGA_param->dimR == 3)
            {
                if(_plyModuleR->getColors().size() == 0)
                    reference_pts = shape_point<double>(_plyModuleR->getVertices(), vel, accl, colorR, p_mass);
                else
                    reference_pts = shape_point<double>(_plyModuleR->getVertices(), vel, accl, _plyModuleR->getColors(), p_mass);
                
                reference_pts.reOrderMissingIndexing(template_pts.getMissingStartIdx(), template_pts.getMissingEndIdx());
                
                NUM_PTS_R = reference_pts.pointSets3D.size();
                
                if(NRGA_param->p_normalize_flag == 1)
                    reference_pts.normalizePoints();
                
                _NRGA->NRGA_initShapePoints_R3D(reference_pts);
                KDTree3D_R = new kdt_nanoFlann<double>(_NRGA->NRGA_getShapePoints_R3D(), reference_pts.pointSets3D.size(), NRGA_param->dimT, 0.001, NRGA_param->k_R);
            }
            else if(NRGA_param->dimR == 4)
            {
                if(_plyModuleR->getColors().size() > 0 && _plyModuleR->getCurvatures().size() > 0)	
                reference_pts = shape_point<double>(_plyModuleR->getVertices(), _plyModuleR->getCurvatures(), vel, accl, _plyModuleR->getColors(), p_mass);
                else
                std::cout << "Vertex Attributes is NOT initialized for REFERENCE" << std::endl;
            }
            
            // save template ply file
            if(NRGA_param->dimT == 3)
            {
                Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_T3D())[0].pos.x,
                    template_pts.pointSets3D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleT->setVertices(T_pts);
            }
            else if(NRGA_param->dimT == 4)
            {
                Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_T4D())[0].pos.x,
                    template_pts.pointSets4D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleT->setVertices(T_pts);
            }
            
            Eigen::MatrixXi Same_color(template_pts.pointSets3D.size(), 3);
            Eigen::VectorXi color_(3); color_ << 0, 255, 0;
            Same_color.setZero();
            Same_color.rowwise() += color_.transpose();
            _plyModuleT->setColors(Same_color);

            QString name = QString(path_Reps.c_str()) + "/" + "Template_Init.ply";  
            bool successT = _plyModuleT->writePLY(name.toStdString(), true, true, false, false, false);
            
            if (successT == true)
                    std::cout <<  "Initial Template file has been read properly" << std::endl;
            
            do
            {
                if(NRGA_param->dimT == 3 && NRGA_param->dimR == 3) 
                {
                    KDTree3D_T = new kdt_nanoFlann<double>(_NRGA->NRGA_getShapePoints_T3D(), template_pts.pointSets3D.size(), NRGA_param->dimT, 0.001, NRGA_param->k_T);
                    
                    _NRGA->NRGA_setCurrIterNo(iter);
                    _NRGA->NRGA_Gravitational_Force3D_ARAP3(_NRGA->NRGA_getShapePoints_T3D(), _NRGA->NRGA_getShapePoints_R3D(), _NRGA->NRGA_getShapePoints_T3DPrev(), KDTree3D_T, KDTree3D_R, template_pts.pointSets3D.size(), reference_pts.pointSets3D.size());
                }
                if(NRGA_param->dimT == 4 && NRGA_param->dimR == 4) 
                {
                    ;
                }
                    
                // save template ply file
// 		    if(NRGA_param->dimT == 3)
// 		    {
// 			Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
// 			    &(_NRGA->NRGA_getShapePoints_T3D())[0].pos.x,
// 			    template_pts.pointSets3D.size(),
// 			    3,
// 			    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
// 			    );
// 			_plyModuleT->setVertices(T_pts);
// 		    }
// 		    else if(NRGA_param->dimT == 4)
// 		    {
// 			Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
// 			    &(_NRGA->NRGA_getShapePoints_T4D())[0].pos.x,
// 			    template_pts.pointSets4D.size(),
// 			    3,
// 			    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
// 			    );
// 			_plyModuleT->setVertices(T_pts);
// 		    }
// 
// 		    name = QString(path_Reps.c_str()) +  + "Template-" + QString::number(iter) + ".ply";	  
// 		    successT = _plyModuleT->writePLY(name.toStdString(), true, true, false, false, false);
            
                std::cout <<  "Reps : " << reps << " ---> Iteration : " << iter++ << std::endl;
                if(iter >= NRGA_param->p_max_iterations)
                    break;
            }
            while(1);
        
            // save Final template ply file-------------------------------------------------------
            if(NRGA_param->dimT == 3)
            {
                Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_T3D())[0].pos.x,
                    template_pts.pointSets3D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleT->setVertices(T_pts);
            }
            else if(NRGA_param->dimT == 4)
            {
                Eigen::MatrixXd T_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_T4D())[0].pos.x,
                    template_pts.pointSets4D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleT->setVertices(T_pts);
            }

            name = QString(path_Reps.c_str()) + "/" + "Template_Final.ply";		  
            successT = _plyModuleT->writePLY(name.toStdString(), true, true, false, false, false);
            
            if (successT == true)
                    std::cout << "Final Template file has been read properly";
            
            //-------------------------------------------------------------------------------------
            
            
            // Saving Final Reference File---------------------------------------------------------
            if(NRGA_param->dimR == 3)
            {
                Eigen::MatrixXd R_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_R3D())[0].pos.x,
                    reference_pts.pointSets3D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point3D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleR->setVertices(R_pts);
            }
            else if(NRGA_param->dimR == 4)
            {
                Eigen::MatrixXd R_pts = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(
                    &(_NRGA->NRGA_getShapePoints_R4D())[0].pos.x,
                    reference_pts.pointSets3D.size(),
                    3,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(offsetof(point4D<double>, m) / sizeof(double) + 1, 1)
                    );
                _plyModuleR->setVertices(R_pts);
            }

            Eigen::MatrixXi Same_colorR(reference_pts.pointSets3D.size(), 3);
            Eigen::VectorXi colorR_(3); color_ << (255), 0, (iter);
            Same_colorR.setZero();
            Same_colorR.rowwise() += colorR_.transpose();
            _plyModuleR->setColors(Same_color);

            name = QString(path_Reps.c_str()) + "/" + "Reference_Final.ply";		  
            bool successR = _plyModuleR->writePLY(name.toStdString(), false, true, false, false, false);
            if (successR == true)
                std::cout << "Reference file has been read properly" << std::endl;
            
            //--------------------------------------------------------------------------------------
            
            
            // calculate mean deviation between Resultant point clouds
            _NRGA->NRGA_calculateMeanDeviation(_NRGA->NRGA_getShapePoints_T3D(), _NRGA->NRGA_getShapePoints_R3D(), NUM_PTS_T);
            // calculate RMSE between the Resultant point clouds
            _NRGA->NRGA_calculateRMSE(_NRGA->NRGA_getShapePoints_T3D(), _NRGA->NRGA_getShapePoints_R3D(), NUM_PTS_T, _NRGA->meanDeviationTR);
                std::string path = NRGA_param->out_registered + "/" + "statistics.txt";
            appendStatistics(path.c_str(), std::to_string(_NRGA->meanDeviationTR) , std::to_string(_NRGA->RMSE_TR));

        }
        while(1);
            
        return 0;
	}
	else
	{
		;
	}
}
#endif



