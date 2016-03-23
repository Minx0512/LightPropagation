#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (PSNR_SRCS
  PSNR.cpp
)




# Create the LPSimulation C++ library.
ADD_LIBRARY (LPPSNR ${PSNR_SRCS})
# TARGET_LINK_LIBRARIES (LPSimulation ${OPENCL_LIBRARIES} ${VTK_LIBRARIES})


