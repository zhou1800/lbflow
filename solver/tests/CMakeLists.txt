find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import numpy"
    RESULT_VARIABLE numpy_status OUTPUT_QUIET ERROR_QUIET)
if(NOT numpy_status EQUAL 0)
    message(FATAL_ERROR "Tests need NumPy: python3 -m pip install -r requirements-test.txt; or configure -DBUILD_TESTING=OFF")
endif()

add_executable(NumericalAudit numerical_audit.cpp)
target_link_libraries(NumericalAudit PRIVATE lbm)

add_test(NAME numerical_moments
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/scripts/ctest_moments.py"
        --driver $<TARGET_FILE:NumericalAudit> --root "${PROJECT_BINARY_DIR}/verification")
add_test(NAME simulator_contracts
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/scripts/ci_smoke.py"
        --solver $<TARGET_FILE:LBFlow> --root "${PROJECT_BINARY_DIR}/verification"
        --mpiexec "${MPIEXEC_EXECUTABLE}" --numproc-flag=${MPIEXEC_NUMPROC_FLAG})
set_tests_properties(numerical_moments simulator_contracts PROPERTIES TIMEOUT 180)
set_tests_properties(simulator_contracts PROPERTIES PROCESSORS 2)
