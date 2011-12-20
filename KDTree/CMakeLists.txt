LINK_DIRECTORIES(/home/doriad/bin/boost/lib)

#####################
INCLUDE_DIRECTORIES(
${VXL_RPL_INCLUDE_DIR}
${VXL_MUL_INCLUDE_DIR}
)

# # Commented due to ModelFile API change
# ADD_EXECUTABLE(TestKDTree TestKDTree.cpp KDTree.cpp
# 
# ${TOOLS}/Tools.cpp
# 
# ${GEOMETRY}/Angles.cpp
# ${GEOMETRY}/Edge.cpp
# ${GEOMETRY}/Triangle.cpp
# ${GEOMETRY}/Color.cpp
# ${GEOMETRY}/Helpers.cpp
# ${GEOMETRY}/Ray.cpp
# ${GEOMETRY}/OrientedPoint.cpp
# ${GEOMETRY}/Geometry.cpp
# 
# ${MODELFILE}/ModelFile.cpp
# ${VXLHELPERS}/VXLHelpers.cpp
# )
# 
# TARGET_LINK_LIBRARIES(TestKDTree ${VXLLibraries} ${VTK_LIBRARIES})
