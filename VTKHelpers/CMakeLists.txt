add_library(libVTKHelpers
VTKHelpers.cpp
Octree.cpp

${VXLHELPERS}/VXLHelpers.cpp

${GEOMETRY}/Geometry.cpp
${GEOMETRY}/Transformation.cpp
${GEOMETRY}/Angles.cpp
${GEOMETRY}/Color.cpp
${GEOMETRY}/Ray.cpp
${GEOMETRY}/OrientedPoint.cpp

${TOOLS}/Tools.cpp
)

TARGET_LINK_LIBRARIES(libVTKHelpers
${VXLLibraries}
${VTK_LIBRARIES}
)

add_subdirectory(Tests)

################################
ADD_EXECUTABLE(CreateVtpLine CreateVtpLine.cpp)
TARGET_LINK_LIBRARIES(CreateVtpLine ${VTK_LIBRARIES})

################################
ADD_EXECUTABLE(CreateVtpLineParts CreateVtpLineParts.cpp)
TARGET_LINK_LIBRARIES(CreateVtpLineParts ${VTK_LIBRARIES})

################################
ADD_EXECUTABLE(ScaleLinePlot ScaleLinePlot.cpp)
