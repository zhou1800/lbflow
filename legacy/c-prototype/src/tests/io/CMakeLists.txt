# libio test executable
add_executable(io-test "io-test.cpp")
target_link_libraries(io-test
    gtest
    libIO
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if ($ENV{N_PROC})
    add_test(NAME io_test COMMAND mpirun -n $ENV{N_PROC} ./io-test)
else()
    add_test(NAME io_test COMMAND io-test)
endif()
