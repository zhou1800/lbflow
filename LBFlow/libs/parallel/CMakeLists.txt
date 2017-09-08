aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
add_library(libParallel
    ${source_files}
)
target_link_libraries(libParallel
    ${MPI_CXX_LIBRARIES}
)
