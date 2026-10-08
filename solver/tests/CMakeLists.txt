#add debugger
add_library(ldebugger debugger)
target_link_libraries(ldebugger grid)
include_directories(./)

#add tests
set(TESTS tgrid tfile talg)

foreach(test ${TESTS})
    # library
    aux_source_directory(${test} ltest_files)
    add_library(l${test} ${ltest_files})
    unset(ltest_files)
    #include
    include_directories(${test}/)
    #exes
    add_executable(${test} ${test}_run.cpp)
    target_link_libraries(${test} l${test} ldebugger)
endforeach(test)
