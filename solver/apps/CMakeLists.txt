#add executables
#set(EXES Connected Histogram LBSim MarkEdge)
set(EXES LBSim)
foreach(exe ${EXES})
    aux_source_directory(${exe} exe_files)
    add_executable(${exe} ${exe_files})
    unset(exe_files)
endforeach(exe)

#linkage
target_link_libraries(LBSim lbm)
#target_link_libraries(Histogram pipt)
#target_link_libraries(MarkEdge pipt)
#target_link_libraries(Connected pipt)
