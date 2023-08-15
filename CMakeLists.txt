CMAKE_POLICY(SET CMP0020 NEW)
Project(NRGA)
cmake_minimum_required(VERSION 2.8.12)
Message(STATUS "-------------------------------")
Message(STATUS "Processing Project ${PROJNAME}:")


#####################################################################################
############### Source headers and all other files for this project #################
#####################################################################################
IF(UNIX)
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
ENDIF(UNIX)

	

	
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
