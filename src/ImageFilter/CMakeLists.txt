#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (ImageFilter_SRCS ImageFilter.cpp)


# set(EXTRA_FILES_SIM )

#FIND_PACKAGE(OpenCL REQUIRED)
#  INCLUDE(${OPENCL_INCLUDE_DIR})


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


# Create the LPSimulation C++ library.
ADD_LIBRARY (LPImageFilter ${ImageFilter_SRCS})
TARGET_LINK_LIBRARIES (LPImageFilter ${VTK_LIBRARIES})


