cmake_minimum_required(VERSION 2.6)

PROJECT(Tools)

set(GEOMETRY ${CMAKE_CURRENT_SOURCE_DIR}/Geometry)
set(VXLHELPERS ${CMAKE_CURRENT_SOURCE_DIR}/VXLHelpers)
set(VTKHELPERS ${CMAKE_CURRENT_SOURCE_DIR}/VTKHelpers)
set(KDTREE ${CMAKE_CURRENT_SOURCE_DIR}/KDTree)
set(MODELFILE ${CMAKE_CURRENT_SOURCE_DIR}/ModelFile)
set(TOOLS ${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(VXL REQUIRED)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(VXLLibraries  vnl vnl_algo vcl vil vbl vgl vgl_algo vul)

add_executable(TestTools TestTools.cpp Tools.cpp )
target_link_libraries(TestTools ${VXLLibraries})

add_executable(TestValidType TestValidType.cpp ValidType.cpp )

add_subdirectory(KDTree)
add_subdirectory(VTKHelpers)
add_subdirectory(VXLHelpers)
