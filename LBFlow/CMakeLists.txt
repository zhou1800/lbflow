cmake_minimum_required(VERSION 2.8)

project(LBFlow C CXX)
SET(CMAKE_CXX_FLAGS "-std=c++11")

# Parallel computing libraries
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# Headers
include_directories(${PROJECT_SOURCE_DIR}/include/)

# Source files
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/)

# Testing
add_subdirectory(${PROJECT_SOURCE_DIR}/tests/)

add_executable(${PROJECT_NAME} "main.c")
