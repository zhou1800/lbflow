# MPI parallel computing libraries
find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})

# The parallel library
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
add_library(libParallel ${source_files})
target_include_directories(libParallel PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_PATH}
)

target_link_libraries(libParallel
    ${MPI_C_LIBRARIES}
)
