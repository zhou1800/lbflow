cmake_minimum_required(VERSION 2.8)

project(LBFlow C CXX)
SET(CMAKE_CXX_FLAGS "-std=c++11")

# Testing
add_subdirectory(${PROJECT_SOURCE_DIR}/tests/)

add_executable(${PROJECT_NAME} "main.c")
