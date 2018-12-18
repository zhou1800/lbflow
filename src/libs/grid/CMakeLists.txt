# The Grid library
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
add_library(libGrid ${source_files})
target_include_directories(libGrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(libGrid
    libParallel
)
