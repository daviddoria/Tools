add_subdirectory(Tests)

ADD_EXECUTABLE(TestVXLHelpers TestVXLHelpers.cpp VXLHelpers.cpp )

target_link_libraries(TestVXLHelpers ${VXLLibraries})
