# libio test executable
add_executable(io-test "io-test.cpp")
target_link_libraries(io-test
    gtest
    libIO
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME io_test COMMAND io-test)
