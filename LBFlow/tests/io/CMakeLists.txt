# libio test executable
add_executable(testio "testio.cpp")
target_link_libraries(testio
    gtest
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Include for libio test only
include_directories(${PROJECT_SOURCE_DIR}/libs/io)

add_test(NAME test_io COMMAND testio)
