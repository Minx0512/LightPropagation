#
# Source files
#
# Here is where you can add the name of your local common classes.
#

SET (IO_SRCS ImportExport.cpp ReadXML.cpp Logs.cpp)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost REQUIRED COMPONENTS regex)
include_directories(${Boost_INCLUDE_DIRS})

# Create the LPIO C++ library.
ADD_LIBRARY (LPIO ${IO_SRCS})
TARGET_LINK_LIBRARIES (LPIO LPImageFilter ${VTK_LIBRARIES} ${Boost_LIBRARIES})

