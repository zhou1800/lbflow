cmake_minimum_required(VERSION 3.4)

project(LBFlow C CXX)
SET(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-lgcov")

# Headers
include_directories(${PROJECT_SOURCE_DIR}/include/)

# Libraries
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/)

# Testing
add_subdirectory(${PROJECT_SOURCE_DIR}/tests/)

# Third party libraries
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/)

add_executable(${PROJECT_NAME} "main.c")
