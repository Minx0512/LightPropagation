
Set(Boost_NO_BOOST_CMAKE ON)


set(LightPropagation_CMAKE_DIR "${LightPropagation_SOURCE_DIR}/CMake")
set(CMAKE_MODULE_PATH ${LightPropagation_CMAKE_DIR} "${LightPropagation_SOURCE_DIR}/config" ${CMAKE_MODULE_PATH})

# Auto-select bitness based on platform
if( NOT BITNESS )
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(BITNESS 64)
    else()
        set(BITNESS 32)
    endif()
endif()

# Select bitness for non-msvc platform. Can be specified as -DBITNESS=32/64 at command-line
if( NOT MSVC )
    set(BITNESS ${BITNESS} CACHE STRING "Specify bitness")
    set_property(CACHE BITNESS PROPERTY STRINGS "64" "32")
endif()
# Unset LIBRARIES, so that corresponding arch specific libs are found when bitness is changed
unset(OPENCL_LIBRARIES CACHE)
unset(GLUT_LIBRARIES CACHE)
unset(GLEW_LIBRARIES CACHE)

if( BITNESS EQUAL 64 )
    set(BITNESS_SUFFIX x86_64)
elseif( BITNESS EQUAL 32 )
    set(BITNESS_SUFFIX x86)
else()
    message( FATAL_ERROR "Bitness specified is invalid" )
endif()

# Set CMAKE_BUILD_TYPE (default = Release)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE Release)
endif()

# Set platform
if( NOT UNIX )
	set(PLATFORM win)
else()
	set(PLATFORM lnx)
endif()

#
# Configure output paths for libraries and executables.
#
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../lib/${PLATFORM}/${BITNESS_SUFFIX}/${CMAKE_BUILD_TYPE} CACHE PATH
    "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../bin/${PLATFORM}/${BITNESS_SUFFIX}/${CMAKE_BUILD_TYPE} CACHE PATH
    "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
############################################################################


#
# Build shared libs ?
#
# Defaults to the same VTK setting.
#

# Standard CMake option for building libraries shared or static by default.
OPTION(BUILD_SHARED_LIBS "Build with shared libraries." ON)
# Copy the CMake option to a setting with VTKMY_ prefix for use in
# our project.  This name is used in vtkmyConfigure.h.in.
SET(LightPropagation_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})



#
# Try to find VTK and include its settings (otherwise complain)
#
IF(NOT VTK_BINARY_DIR)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
ENDIF(NOT VTK_BINARY_DIR)


  FIND_PACKAGE(OpenCL REQUIRED)
  INCLUDE(${OPENCL_INCLUDE_DIR})



FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )



############################################################################
# Find GL, Glew and Glut libraries

find_path( GL_INCLUDE_DIRS
            NAMES GL/glew.h GL/glut.h
            HINTS ../../../../include/ $ENV{AMDAPPSDKROOT}/include/
)
mark_as_advanced(GL_INCLUDE_DIRS)

find_library( GLUT_LIBRARIES
	NAMES glut glut${BITNESS} freeglut
	HINTS ../../../../lib/ $ENV{AMDAPPSDKROOT}/lib/
	PATH_SUFFIXES ${PLATFORM}${BITNESS} ${BITNESS_SUFFIX}
)
mark_as_advanced( GLUT_LIBRARIES )

find_library( GLEW_LIBRARIES
	NAMES glew GLEW glew${BITNESS} GLEW${BITNESS}
	HINTS ../../../../lib/ $ENV{AMDAPPSDKROOT}/lib/
	PATH_SUFFIXES ${PLATFORM}${BITNESS} ${BITNESS_SUFFIX}
)
mark_as_advanced( GLEW_LIBRARIES )

if( GL_INCLUDE_DIRS STREQUAL "" OR GLUT_LIBRARIES STREQUAL "" OR GLEW_LIBRARIES STREQUAL "")
	message( FATAL_ERROR "Could not locate Glut, Glew include & libs" )
endif( )



set( COMPILER_FLAGS " " )
set( LINKER_FLAGS " " )
set( ADDITIONAL_LIBRARIES "" )

file(GLOB INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.h" )
include_directories( ${OPENCL_INCLUDE_DIR} ${GL_INCLUDE_DIRS} ../../../../include/SDKUtil $ENV{AMDAPPSDKROOT}/include/SDKUtil )

SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

# gcc/g++ specific compile options
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set( COMPILER_FLAGS "${COMPILER_FLAGS} -msse2 " )
    
    # Note: "rt" is not present on mingw
    if( UNIX )
		if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
			set( COMPILER_FLAGS " -g " )
		endif( )
        set( ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} "rt" )
    endif( )
    
    if( BITNESS EQUAL 32 )
        set( COMPILER_FLAGS "${COMPILER_FLAGS} -m32 " )
        set( LINKER_FLAGS "${LINKER_FLAGS} -m32 " )
    else( )
        set( COMPILER_FLAGS "${COMPILER_FLAGS} -m64 " )
        set( LINKER_FLAGS "${LINKER_FLAGS} -m64 " )
    endif( )
    
    # Link GL library. cygwin, mingw names GL libs as libopengl32, linux names it as libGL
    if( MINGW OR CYGWIN )
        set( ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} "opengl32" "glu32" )
    else( )     # Native Linux OS
        set( ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} "GL" )
    endif( )
    
    set( COMPILER_FLAGS "${COMPILER_FLAGS} ${EXTRA_COMPILER_FLAGS_GXX} " )
    set( LINKER_FLAGS "${LINKER_FLAGS} ${EXTRA_LINKER_FLAGS_GXX} " )
    set( ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${EXTRA_LIBRARIES_GXX} )
elseif( MSVC )
    # Samples can specify additional libs/flags using EXTRA* defines
	add_definitions( "/W3 /D_CRT_SECURE_NO_WARNINGS /wd4005 /wd4996 /nologo" )

    set( COMPILER_FLAGS "${COMPILER_FLAGS} ${EXTRA_COMPILER_FLAGS_MSVC} " )
    set( LINKER_FLAGS "${LINKER_FLAGS} ${EXTRA_LINKER_FLAGS_MSVC} /SAFESEH:NO " )
    set( ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${EXTRA_LIBRARIES_MSVC} )
endif( )

#--------------
#....
add_custom_target(copy_patterns ALL
  COMMAND ${CMAKE_COMMAND} -DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DDSTDIR=${CMAKE_CURRENT_BINARY_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/copy_patterns.cmake
  COMMENT "Copying pattern files to build tree")
#....
#--------------------


# Copy extra files to binary directory
#foreach( extra_file ${EXTRA_FILES} )
#    add_custom_command(
#        TARGET LightPropagation POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        ${CMAKE_CURRENT_SOURCE_DIR}/${extra_file}  ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}
#		COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        ${CMAKE_CURRENT_SOURCE_DIR}/${extra_file}  ./
#        )
#endforeach( extra_file )




