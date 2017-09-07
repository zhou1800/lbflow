# io test executable
add_executable(testio "testio.cpp")
target_link_libraries(testio
    gtest
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME test_io COMMAND testio)
