
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


ADD_EXECUTABLE(LightPropagation main.cpp)
TARGET_LINK_LIBRARIES(LightPropagation LPIO LPImageFilter LPModel LPSimulation LPAnalysis LPInterface ${VTK_LIBRARIES})


# set_target_properties( ${SAMPLE_NAME} PROPERTIES
#                      COMPILE_FLAGS ${COMPILER_FLAGS}
#                     LINK_FLAGS ${LINKER_FLAGS}
#                 )
#target_link_libraries( ${SAMPLE_NAME} ${VTK_LIBRARIES} ${OPENCL_LIBRARIES}
#                        ${GLUT_LIBRARIES} ${GLEW_LIBRARIES} ${ADDITIONAL_LIBRARIES} )

