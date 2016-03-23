cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)


set(SAMPLE_NAME LightPropagation)
#
# Set the project name.
#

PROJECT (LightPropagation)

#
# Load CMake commands that you probably should not modify.
#

INCLUDE (${LightPropagation_SOURCE_DIR}/CMakeOptions.cmake)


INCLUDE_DIRECTORIES(${LightPropagation_BINARY_DIR})

#
# You can put your include path(s) here
#

INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/IO)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/Model)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/Simulation)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/PSNR)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/Analysis)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/Interface)
INCLUDE_DIRECTORIES(${LightPropagation_SOURCE_DIR}/ImageFilter)



# set( EXTRA_FILES ${EXTRA_FILES_SIM} ${EXTRA_FILES_MODEL})

#
# Build main programm  too ?
#

OPTION(BUILD_MAIN "Build Programm." ON)
IF (BUILD_MAIN)
  ADD_SUBDIRECTORY(main)
ENDIF (BUILD_MAIN)



#
# Here is where you can list the sub-directories holding your local
# classes. Sorting classes by 'package' type like VTK does (Common,
# Rendering, Filtering, Imaging, IO, etc.) is a good thing and prevents
# numerous dependencies problems.
#

ADD_SUBDIRECTORY(IO)
ADD_SUBDIRECTORY(Model)
ADD_SUBDIRECTORY(Simulation)
ADD_SUBDIRECTORY(PSNR)
ADD_SUBDIRECTORY(Analysis)
ADD_SUBDIRECTORY(Interface)
ADD_SUBDIRECTORY(ImageFilter)

# set(EXTRA_FILES_MODEL layermodel.xml)
# set(EXTRA_FILES ${EXTRA_FILES_MODEL})

# set(EXTRA_FILES ${EXTRA_FILES_MODEL} ${EXTRA_FILES_SIM})


