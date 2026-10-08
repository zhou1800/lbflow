#add .h files to CMake project
FILE(GLOB_RECURSE HFiles "*.h")
add_custom_target(headers SOURCES ${HFiles})
