#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Interface_SRCS
  CLI.cpp
)




# Create the LPSimulation C++ library.
ADD_LIBRARY (LPInterface ${Interface_SRCS})
# TARGET_LINK_LIBRARIES (LPSimulation ${OPENCL_LIBRARIES} ${VTK_LIBRARIES})


