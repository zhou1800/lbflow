#add libraries
set(MODS algorithm file grid helper lbm parallel pipt)
foreach(mod ${MODS})
    aux_source_directory(${mod} mod_files)
    add_library(${mod} ${mod_files})
    unset(mod_files)
endforeach(mod)

#linkage
target_link_libraries(parallel ${MPI_CXX_LIBRARIES})

target_link_libraries(grid helper)
target_link_libraries(grid algorithm)
target_link_libraries(file parallel)
target_link_libraries(helper file)
target_link_libraries(lbm grid file)
target_link_libraries(pipt grid file)
