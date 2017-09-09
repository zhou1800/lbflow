# libparallel test executable
add_executable(testparallel "testparallel.cpp")
target_link_libraries(testparallel
    gtest
    libParallel
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME test_parallel COMMAND testparallel)
