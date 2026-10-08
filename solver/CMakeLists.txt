cmake_minimum_required(VERSION 2.8)
project(DRLab CXX)
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS "-std=c++11")

find_package(MPI REQUIRED)

# set directory
set(SRC_DIR src/)
set(APP_DIR apps/)
set(TEST_DIR tests/)
set(INC_DIR ${CMAKE_SOURCE_DIR}/inc/)

if(MSVC)
    set(IS_MSVC 1)
else(MSVC)
    set(IS_MSVC 0)
endif(MSVC)

#add configure
configure_file(
"${INC_DIR}core/cfg.h.in"
"${INC_DIR}core/cfg.h"
)

#add includes
include_directories(${MPI_CXX_INCLUDE_PATH})
set(INCS algorithm core file grid helper lbm parallel pipt)
foreach(inc ${INCS})
    include_directories(${INC_DIR}${inc}/)
endforeach(inc)

#add sub directory
add_subdirectory(${SRC_DIR})
add_subdirectory(${APP_DIR})
add_subdirectory(${TEST_DIR})
add_subdirectory(${INC_DIR})
