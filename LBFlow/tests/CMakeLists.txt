enable_testing()

# Google test
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/googletest/)

# All tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io/)
