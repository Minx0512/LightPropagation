#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Simulation_SRCS
  CalculatePropagation.cpp
)


# set(EXTRA_FILES_SIM MCLightPropagation_CL.cl vectorAddKernel.cl)

FIND_PACKAGE(OpenCL REQUIRED)
  INCLUDE(${OPENCL_INCLUDE_DIR})


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


# Create the LPSimulation C++ library.
ADD_LIBRARY (LPSimulation ${Simulation_SRCS})
TARGET_LINK_LIBRARIES (LPSimulation LPPSNR ${OPENCL_LIBRARIES} ${VTK_LIBRARIES})


