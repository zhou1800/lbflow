# grid test executable
add_executable(grid-test "grid-test.cpp")
target_link_libraries(grid-test
    gtest
    libGrid
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if ($ENV{N_PROC})
    add_test(NAME grid_test COMMAND mpirun -n $ENV{N_PROC} ./grid-test)
else()
    add_test(NAME grid_test COMMAND grid-test)
endif()
