ENABLE_TESTING()

ADD_EXECUTABLE(TestRotation TestRotation.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ../Angles.cpp ${PROJECTS}/Geometry/Transformation.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestRotation ${Libraries})
ADD_TEST(TestRotation TestRotation)

ADD_EXECUTABLE(TestFitPlane TestFitPlane.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestFitPlane ${Libraries} libGeometry)
ADD_TEST(TestFitPlane TestFitPlane)

ADD_EXECUTABLE(TestAligningRotation TestAligningRotation.cpp)
TARGET_LINK_LIBRARIES(TestAligningRotation ${Libraries})
ADD_TEST(TestAligningRotation TestAligningRotation)

ADD_EXECUTABLE(TestBestPlane TestBestPlane.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestBestPlane ${Libraries} libGeometry)
ADD_TEST(TestBestPlane TestBestPlane)

ADD_EXECUTABLE(TestRayIsSameDirection TestRayIsSameDirection.cpp ${PROJECTS}/Geometry/Ray.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestRayIsSameDirection ${Libraries})
ADD_TEST(TestRayIsSameDirection TestRayIsSameDirection)

ADD_EXECUTABLE(TestCorrespondenceError TestCorrespondenceError.cpp)
TARGET_LINK_LIBRARIES(TestCorrespondenceError ${Libraries})
ADD_TEST(TestCorrespondenceError TestCorrespondenceError)

ADD_EXECUTABLE(TestIsRightHanded TestIsRightHanded.cpp ${PROJECTS}/Geometry/Angles.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestIsRightHanded ${Libraries})
ADD_TEST(TestIsRightHanded TestIsRightHanded)

ADD_EXECUTABLE(TestLineSegmentRotation TestLineSegmentRotation.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestLineSegmentRotation ${Libraries} libGeometry)
ADD_TEST(TestLineSegmentRotation TestLineSegmentRotation)

ADD_EXECUTABLE(TestTransformation TestTransformation.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Geometry/Transformation.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestTransformation ${Libraries})
ADD_TEST(TestTransformation TestTransformation)

ADD_EXECUTABLE(TestAxisAngle TestAxisAngle.cpp)
TARGET_LINK_LIBRARIES(TestAxisAngle ${Libraries})
ADD_TEST(TestAxisAngle TestAxisAngle)

ADD_EXECUTABLE(TestIntersectRayBox TestIntersectRayBox.cpp ${PROJECTS}/Geometry/Ray.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp ${PROJECTS}/Tools/Tools.cpp)
TARGET_LINK_LIBRARIES(TestIntersectRayBox ${Libraries})
ADD_TEST(TestIntersectRayBox TestIntersectRayBox)

ADD_EXECUTABLE(TestSphereRect TestSphereRect.cpp  ${PROJECTS}/Geometry/Angles.cpp ${PROJECTS}/VXLHelpers/VXLHelpers.cpp)
TARGET_LINK_LIBRARIES(TestSphereRect ${Libraries})
ADD_TEST(TestSphereRect TestSphereRect)
