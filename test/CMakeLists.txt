
function(setup_test_dir TESTDIR)
  # make test directory within build tree (will be used for test output/current working directory)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  # link in input directory because of fixed input file hierarchy 
  file(CREATE_LINK ${PROJECT_SOURCE_DIR}/input ${PROJECT_BINARY_DIR}/test/${TESTDIR}/input SYMBOLIC)
endfunction()

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

setup_test_serial(Heat-1 MethodHeat/Heat-1  input/Heat/method_heat-1.in)
setup_test_parallel(Heat-8 MethodHeat/Heat-8  input/Heat/method_heat-8.in)

