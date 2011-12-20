ENABLE_TESTING()

ADD_EXECUTABLE(TestEigen3x3 TestEigen3x3.cpp ../VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestEigen3x3 ${VXLLibraries} )
ADD_TEST(TestEigen3x3 TestEigen3x3)

ADD_EXECUTABLE(TestEigenBig TestEigenBig.cpp ../VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestEigenBig ${VXLLibraries})
ADD_TEST(TestEigenBig TestEigenBig)

ADD_EXECUTABLE(TestVectorize TestVectorize.cpp ../VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestVectorize ${VXLLibraries})
ADD_TEST(TestVectorize TestVectorize)

ADD_EXECUTABLE(TestReshape TestReshape.cpp ../VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestReshape ${VXLLibraries})
ADD_TEST(TestReshape TestReshape)

ADD_EXECUTABLE(TestMatrixPowers TestMatrixPowers.cpp ../VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestMatrixPowers ${VXLLibraries})
ADD_TEST(TestMatrixPowers TestMatrixPowers)