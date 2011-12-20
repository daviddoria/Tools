add_library(libGeometry
Angles.cpp
Color.cpp
Edge.cpp
Ray.cpp
Triangle.cpp
Transformation.cpp
LineSegment.cpp
OrientedPoint.cpp
Geometry.cpp
Helpers.cpp
)

ADD_EXECUTABLE(Test Test.cpp
${GeometrySources}
${PROJECTS}/VXLHelpers/VXLHelpers.cpp
${PROJECTS}/Tools/Tools.cpp
#${PROJECTS}/KDTree/KDTree.cpp
#${PROJECTS}/Octree/Octree.cpp
#${PROJECTS}/ModelFile/ModelFile.cpp
)
target_link_libraries(Test ${VTK_LIBRARIES} ${VXLLibraries})

add_subdirectory(Tests)



