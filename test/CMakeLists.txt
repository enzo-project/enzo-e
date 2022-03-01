
function(setup_test_dir TESTDIR)
  # make test directory within build tree (will be used for test output/current working directory)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  # link in input directory because of fixed input file hierarchy 
  file(CREATE_LINK ${PROJECT_SOURCE_DIR}/input ${PROJECT_BINARY_DIR}/test/${TESTDIR}/input SYMBOLIC)
endfunction()

# Test wrapper functions calling enzo-e directly
function(setup_test_serial TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND enzo-e ${INFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial" )
endfunction()

function(setup_test_parallel TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${PARALLEL_LAUNCHER} ${PARALLEL_LAUNCHER_NPROC_ARG} ${PARALLEL_LAUNCHER_NPROC} ${PROJECT_BINARY_DIR}/bin/enzo-e ${INFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel" )
endfunction()

# Test wrapper functions calling a script (e.g., to include output analysis)
function(setup_test_parallel_script TESTNAME TESTDIR TESTCMD TESTSCRIPT)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${TESTCMD} ${TESTSCRIPT} --launch_cmd="${PARALLEL_LAUNCHER} ${PARALLEL_LAUNCHER_NPROC_ARG} ${PARALLEL_LAUNCHER_NPROC} ${PROJECT_BINARY_DIR}/bin/enzo-e"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel" )
endfunction()


setup_test_serial(Heat-1 MethodHeat/Heat-1  input/Heat/method_heat-1.in)
setup_test_parallel(Heat-8 MethodHeat/Heat-8  input/Heat/method_heat-8.in)

setup_test_parallel_script(dual_energy_shock_tube vlct python "input/vlct/run_dual_energy_shock_tube_test.py")

# Convert markdown file to html file for more flexible viewing
configure_file(${CMAKE_SOURCE_DIR}/test/TEST_RESULTS.md ${PROJECT_BINARY_DIR}/TEST_RESULTS.md COPYONLY)
add_custom_target(process_test_results
                  pandoc -s -o TEST_RESULTS.html --metadata title="Enzo-E Test Results" TEST_RESULTS.md
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

