# The IO library
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} source_files)
add_library(libIO ${source_files})
target_include_directories(libIO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(libIO
    cJSON
    libParallel
)
