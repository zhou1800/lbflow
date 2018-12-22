# libparallel test executable
add_executable(parallel-test "parallel-test.cpp")
target_link_libraries(parallel-test
    gtest
    libParallel
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if ($ENV{N_PROC})
    add_test(NAME parallel_test COMMAND mpirun -n $ENV{N_PROC} ./parallel-test)
else()
    add_test(NAME parallel_test COMMAND parallel-test)
endif()
