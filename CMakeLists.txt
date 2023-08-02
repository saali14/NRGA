CMAKE_POLICY(SET CMP0020 NEW)
Project(NRGA)
cmake_minimum_required(VERSION 2.8.12)
Message(STATUS "-------------------------------")
Message(STATUS "Processing Project ${PROJNAME}:")


#####################################################################################
############### Source headers and all other files for this project #################
#####################################################################################
IF(WIN32 OR UNIX)

FILE(GLOB SOURCE_FILES 			${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                                        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
FILE(GLOB HEADER_FILES 			${CMAKE_CURRENT_SOURCE_DIR}/headers/*.hpp
                                        ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h )
FILE(GLOB HELPER_SOURCE_FILES           ${CMAKE_CURRENT_SOURCE_DIR}/helpers/*.cpp
                                        ${CMAKE_CURRENT_SOURCE_DIR}/helpers/*.hpp
                                        ${CMAKE_CURRENT_SOURCE_DIR}/helpers/*.inl
                                        ${CMAKE_CURRENT_SOURCE_DIR}/helpers/*.h
                                        ${CMAKE_CURRENT_SOURCE_DIR}/helpers/*.c)

FILE(GLOB GLSL_FILES 			${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.vert
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.frag
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.geom
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.cont
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.eval
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.comp
                                        ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.glsl)

SOURCE_GROUP(SHADER FILES 		${GLSL_FILES})
SOURCE_GROUP(HELPER FILES		${HELPER_SOURCE_FILES})
ENDIF(WIN32 OR UNIX)


#######################################  WINDOWS OPERATING SYSTEM ##########################################################
############################################################################################################################
IF(WIN32)

#####################################################################################
############### Defining the external library related bariables #####################
#####################################################################################
#SET(CMAKE_CXX_FLAGS 			"${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
SET(CMAKE_PREFIX_PATH 			"C:\\Qt\\Qt5.5.1\\5.5\\msvc2013_64\\")
SET(EIGEN_DIR                           ${CMAKE_CURRENT_SOURCE_DIR}/libs/Eigen_source)
SET(GLEW_INCLUDE_DIR                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/glew2_0/include)
SET(GLEW_LIB_32_DIR                     ${CMAKE_CURRENT_SOURCE_DIR}/libs/glew2_0/lib/Release/Win32)
SET(GLEW_LIB_64_DIR                     ${CMAKE_CURRENT_SOURCE_DIR}/libs/glew2_0/lib/Release/x64)
#SET(KDTREE_ANN_INCLUDE_DIR              ${CMAKE_CURRENT_SOURCE_DIR}/libs/ann/include)
#SET(KDTREE_ANN_LIB_DIR                  ${CMAKE_CURRENT_SOURCE_DIR}/libs/ann/lib)
#SET(KDTREE_ANN_BIN_DIR                  ${CMAKE_CURRENT_SOURCE_DIR}/libs/ann/bin)
SET(FLANN_INCLUDE_DIRS			"C:/Program Files/flann/include")
SET(FLANN_LIBRARY_DIRS			"C:/Program Files/flann/lib")
SET(LIBCOFIGPP_INCLUDE_DIR		"C:/Program Files (x86)/GnuWin32/include" CACHE PATH  "Libconfig Include Path")
SET(LIBCOFIGPP_LIB_DIR			"C:/Program Files (x86)/GnuWin32/lib" CACHE PATH "Libconfig Library Path")
SET(LIBCOFIGPP_BIN_DIR			"C:/Program Files (x86)/GnuWin32/bin" CACHE PATH "Libconfig Bin Path")
SET(NANOFLANN_INCLUDE_DIRS		${CMAKE_CURRENT_SOURCE_DIR}/libs/nanoflann/include)
SET(PLY_LIB                             ${CMAKE_CURRENT_SOURCE_DIR}/libs)
SET(CMAKE_INCLUDE_PATH 			${CMAKE_INCLUDE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0/")
SET(CMAKE_LIBRARY_PATH 			${CMAKE_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0/stage/x64/lib")
SET(BOOST_INCLUDE_DIR			${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0)
SET(BOOST_ROOT				${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0)
SET(BOOST_LIB_DIR			${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0/stage/x64/lib)



#####################################################################################
################ Finding the Libraries packages from modules ########################
#####################################################################################
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(OpenGL)
FIND_PACKAGE(Qt5Widgets)
FIND_PACKAGE(Qt5Core)
FIND_PACKAGE(Qt5Gui)
FIND_PACKAGE(Qt5OpenGL)



#####################################################################################
########### Adding the definitions of the packages found by module system ###########
#####################################################################################
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(${Qt5Widget_DEFINITIONS})


#SET(CUDA_NVCC_FLAGS "-arch=sm_20" CACHE STRING "nvcc flags" FORCE)
#SET(CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)
#CUDA_ADD_EXECUTABLE(cuda_NRGA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
#CUDA_ADD_LIBRARY(cuda_NRGA_lib ${CUDA_FILES})
#ADD_DEPENDENCIES(cuda_NRGA cuda_NRGA_lib)
#TARGET_LINK_LIBRARIES(cuda_NRGA cuda_NRGA_lib)


#####################################################################################
###### Let compiler know which directories to search for " Including " headers ######
#####################################################################################
INCLUDE_DIRECTORIES(
                                        ${CMAKE_CURRENT_SOURCE_DIR}/headers
                                        #${KDTREE_ANN_INCLUDE_DIR}
                                        ${GLEW_INCLUDE_DIR}
                                        #${SOIL_INCLUDE_DIR}
                                        ${LIBCOFIGPP_INCLUDE_DIR}
                                        ${AUTODIFF_INCLUDE_DIR}
                                        ${EIGEN_DIR}
                                        ${Boost_INCLUDE_DIR}
                                        ${BOOST_INCLUDE_DIR}
                                        ${FLANN_INCLUDE_DIRS}
                                        ${NANOFLANN_INCLUDE_DIRS}
)


#####################################################################################
###### Let compiler know which directories to link against for the libraries ########
#####################################################################################
LINK_DIRECTORIES(
                                        ${PLY_LIB}
                                        ${GLEW_LIB_64_DIR}
                                        #${KDTREE_ANN_LIB_DIR}
                                        #${KDTREE_ANN_BIN_DIR}
                                        ${LIBCOFIGPP_LIB_DIR}
                                        ${LIBCOFIGPP_BIN_DIR}
					${Qt5Core_QTMAIN_LIBRARIES}
                                        ${Boost_LIBRARY_DIRS}
                                        ${BOOST_LIB_DIR}
                                        ${FLANN_LIBRARY_DIRS}
)


#####################################################################################
############################### Executable Linkage ##################################
#####################################################################################
ADD_EXECUTABLE(${PROJECT_NAME} main.cpp ${HEADER_FILES} 
					${SOURCE_FILES} 
					${HELPER_SOURCE_FILES} 
#					${CUDA_FILES} 
					${GLSL_FILES}
)
qt5_use_modules(${PROJECT_NAME} Widgets Core OpenGL)

target_link_libraries(${PROJECT_NAME}	 
					${CUDA_LIBRARIES}
					${LIBCONFIGPP_LIBRARY} 
					${Boost_LIBRARIES} 
					${QT_LIBRARIES} 
					${OPENGL_LIBRARIES}
					#SOIL 
					glew32s 
					#ANN 
#					cuda_NRGA_lib 
					flann
)
ENDIF(WIN32)

	

	
#######################################  LINUX OPERATING SYSTEM ##########################################################
##########################################################################################################################
IF(UNIX)

# guard for adding different compilation type
if(CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE(STATUS "GCC detected, adding compile flags")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -fPIC")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -fPIC")

#####################################################################################
############### Defining the external library related bariables #####################
#####################################################################################
SET(EIGEN_DIR                           ${CMAKE_CURRENT_SOURCE_DIR}/libs/Eigen_source)
SET(GLEW_INCLUDE_DIR                    ${CMAKE_CURRENT_SOURCE_DIR}/libs/glew2_0/include)
SET(GLEW_LIB_DIR                     	${CMAKE_CURRENT_SOURCE_DIR}/libs/glew2_0/lib/)
SET(FLANN_INCLUDE_DIRS			${FLANN_INCLUDE_DIR})
SET(FLANN_LIBRARY_DIRS			${FLANN_LIBRARY})
SET(NANOFLANN_INCLUDE_DIRS		${CMAKE_CURRENT_SOURCE_DIR}/libs/nanoflann/include)
SET(CMAKE_INCLUDE_PATH 			${CMAKE_INCLUDE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0/")
SET(CMAKE_LIBRARY_PATH 			${CMAKE_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0/stage/x64/lib")
SET(BOOST_ROOT				${CMAKE_CURRENT_SOURCE_DIR}/libs/boost_1_58_0)
SET(CMAKE_AUTOMOC 			ON)
SET(CMAKE_AUTOUIC 			ON)
SET(CMAKE_INCLUDE_CURRENT_DIR 		ON)
set(QT_USE_QTOPENGL 			TRUE)

#FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)
FIND_PACKAGE(Qt5Widgets REQUIRED)
FIND_PACKAGE(Qt5Core REQUIRED)
FIND_PACKAGE(Qt5Gui REQUIRED)
FIND_PACKAGE(Qt5OpenGL REQUIRED)

IF(NOT OPENGL_FOUND)
MESSAGE(ERROR " OPENGL not found!")
ENDIF(NOT OPENGL_FOUND)
IF(NOT GLUT_FOUND)
MESSAGE(ERROR " GLUT not found!")
ENDIF(NOT GLUT_FOUND)

IF(Qt5_POSITION_INDEPENDENT_CODE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

FIND_PATH(LIBCONFIG_INCLUDE_DIR libconfig.h /usr/local/include /usr/include ~/libs/include ~/include)
FIND_PATH(LIBCONFIGPP_INCLUDE_DIR libconfig.h++ /usr/local/include /usr/include ~/libs/include ~/include)
FIND_LIBRARY(LIBCONFIGPP_LIBRARY config++ /usr/local/lib /usr/lib ~/libs/lib)
FIND_LIBRARY(LIBCONFIG_LIBRARY config /usr/local/lib /usr/lib ~/libs/lib)
IF(LIBCONFIGPP_INCLUDE_DIR AND LIBCONFIGPP_LIBRARY)
SET(LIBCONFIGPP_FOUND TRUE)
ENDIF(LIBCONFIGPP_INCLUDE_DIR AND LIBCONFIGPP_LIBRARY)
IF(LIBCONFIG_INCLUDE_DIR AND LIBCONFIG_LIBRARY)
SET(LIBCONFIG_FOUND TRUE)
ENDIF(LIBCONFIG_INCLUDE_DIR AND LIBCONFIG_LIBRARY)
	
ADD_DEFINITIONS(${OpenGL_DEFINITIONS})
ADD_DEFINITIONS(${GLUT_DEFINITIONS})
ADD_DEFINITIONS(${Qt5Widget_DEFINITIONS})


#SET(CUDA_NVCC_FLAGS "-arch=sm_20;-std=c++11;-O2;-DVERBOSE" CACHE STRING "nvcc flags" FORCE)
#SET(CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)
#CUDA_ADD_EXECUTABLE(cuda_NRGA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
#CUDA_ADD_LIBRARY(cuda_NRGA_lib ${CUDA_FILES})
#ADD_DEPENDENCIES(cuda_NRGA cuda_NRGA_lib)
#TARGET_LINK_LIBRARIES(cuda_NRGA cuda_NRGA_lib)

INCLUDE_DIRECTORIES(            
					${CMAKE_CURRENT_SOURCE_DIR}/headers
					${Boost_INCLUDE_DIRS}
					${Qt5Widgets_INCLUDE_DIRS}
					${GLEW_INCLUDE_DIR}
					${OPENGL_INCLUDE_DIRS}
					${GLUT_INCLUDE_DIRS}
					${LIBCOFIGPP_INCLUDE_DIR}
					${EIGEN_DIR}
					${FLANN_INCLUDE_DIRS}
					${NANOFLANN_INCLUDE_DIRS}
					)

LINK_DIRECTORIES(		
					${BOOST_LIB}
					${EIGEN_LIB}
					${GLUT_LIBRARY_DIRS}
					${OpenGL_LIBRARY_DIRS}
					${GLEW_LIB_DIR}
					${LIBCONFIGPP_LIBRARY}
					${Qt5Core_QTMAIN_LIBRARIES}
					${FLANN_LIBRARY_DIRS}
					)

ADD_EXECUTABLE( ${PROJECT_NAME}	main.cpp ${NRGA_wrap_src} 
					 ${SOURCE_FILES} 
					 ${HEADER_FILES} 
					 ${HELPER_SOURCE_FILES} 
					 ${GLSL_FILES}
		)
qt5_use_modules(${PROJECT_NAME} Core Widgets OpenGL)
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
					${LIBCONFIGPP_LIBRARY}
					${OPENGL_LIBRARIES}
					${GLUT_LIBRARY}
					${QT_LIBRARIES}
					Qt5::Widgets
					libboost_iostreams.so
					libboost_system.so
					libboost_filesystem.so
					-lGL
					-lGLEW
					-lm
					-lpthread
					)
ENDIF(UNIX)
