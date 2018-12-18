enable_testing()

# All tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/parallel/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/grid/)

# configure file
configure_file(
"${CMAKE_CURRENT_SOURCE_DIR}/testcfg.h.in"
"${CMAKE_CURRENT_SOURCE_DIR}/testcfg.h"
)
