cmake_minimum_required(VERSION 3.5)

project(LBFlow C CXX)
SET(CMAKE_CXX_FLAGS "-std=c++11")

# Headers
include_directories(${PROJECT_SOURCE_DIR}/include/)

# Libraries
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/)

# Testing
add_subdirectory(${PROJECT_SOURCE_DIR}/tests/)

# Third party libraries
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/)

add_executable(${PROJECT_NAME} "main.c")
