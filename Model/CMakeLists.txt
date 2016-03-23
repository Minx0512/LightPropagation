#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Model_SRCS
	Emitter.cpp
	createModel.cpp
	VoxToMesh.cpp
)


# set(EXTRA_FILES_MODEL layermodel.xml)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Create the LPModel C++ library.
ADD_LIBRARY (LPModel ${Model_SRCS})

TARGET_LINK_LIBRARIES (LPModel ${VTK_LIBRARIES})

