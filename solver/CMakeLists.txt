# The repository-root project owns solver configuration.
if(MSVC)
    set(IS_MSVC 1)
else()
    set(IS_MSVC 0)
endif()
configure_file(inc/core/cfg.h.in "${CMAKE_CURRENT_BINARY_DIR}/include/cfg.h")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
foreach(module algorithm core file grid helper lbm parallel pipt)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc/${module}")
endforeach()

add_subdirectory(src)
add_subdirectory(apps)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
