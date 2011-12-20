ENABLE_TESTING()

ADD_EXECUTABLE(TestBestFitPlane TestBestFitPlane.cpp)
TARGET_LINK_LIBRARIES(TestBestFitPlane libVTKHelpers ${VTK_LIBRARIES} ${VXLLibraries})
ADD_TEST(TestBestFitPlane TestBestFitPlane)

#ADD_EXECUTABLE(TestFitPlane TestFitPlane.cpp)
#ADD_TEST(TestFitPlane TestFitPlane)