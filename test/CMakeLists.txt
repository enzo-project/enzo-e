
function(setup_test_dir TESTDIR)
  # make test directory within build tree (will be used for test output/current working directory)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  # link in input directory because of fixed input file hierarchy 
  file(CREATE_LINK ${PROJECT_SOURCE_DIR}/input ${PROJECT_BINARY_DIR}/test/${TESTDIR}/input SYMBOLIC)
endfunction()

# Test wrapper function calling unit tests (separate binaries) directly
function(setup_test_unit TESTNAME TESTDIR TESTBIN)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${TESTBIN}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;unit" )
endfunction()

# Test wrapper functions calling enzo-e directly
function(setup_test_serial TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND enzo-e ${INFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;enzo" )
endfunction()

function(setup_test_parallel TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${PARALLEL_LAUNCHER} ${PARALLEL_LAUNCHER_NPROC_ARG} ${PARALLEL_LAUNCHER_NPROC} ${PROJECT_BINARY_DIR}/bin/enzo-e ${INFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel;enzo" )
endfunction()

# Serial (single compute element) test wrapper functions calling a script (e.g., to include output analysis)
# All optional argment (i.e., the ones not being named) are automatically passed to the command being called.
function(setup_test_serial_python TESTNAME TESTDIR TESTSCRIPT)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${TESTSCRIPT} --launch_cmd=${PROJECT_BINARY_DIR}/bin/enzo-e ${ARGN}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;enzo;regression" )
endfunction()

# Parallel (using multiple compute elements) test wrapper functions calling a script (e.g., to include output analysis)
# All optional argment (i.e., the ones not being named) are automatically passed to the command being called.
function(setup_test_parallel_python TESTNAME TESTDIR TESTSCRIPT)
  # construct a list of the various parts of the launch_cmd argument (this is
  # necessary since PARALLEL_LAUNCHER_NPROC_ARG may itself be a list)
  set(LAUNCH_CMD_LIST
    ${PARALLEL_LAUNCHER}
    ${PARALLEL_LAUNCHER_NPROC_ARG} ${PARALLEL_LAUNCHER_NPROC}
    ${PROJECT_BINARY_DIR}/bin/enzo-e)
  # now join the elements of LAUNCH_CMD_LIST into a single string so that the
  # shell interprets the entire string as a single command line argument
  LIST(JOIN LAUNCH_CMD_LIST " " LAUNCH_CMD)

  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${TESTSCRIPT} --launch_cmd ${LAUNCH_CMD} ${ARGN}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel;enzo;regression" )
endfunction()

############################### UNIT TESTS ####################################
# The following tests are self-contained binary of individual components

#setup_test_unit(Array ArrayComponent/Array test_cello_array)
setup_test_unit(Memory MemoryComponent/Memory test_memory)
setup_test_unit(Monitor MonitorComponent/Monitor test_monitor)
#setup_test_unit( Component/ test_)

############################### ENZO-E TESTS ##################################
# The following tests will call the enzo-e binary in one way or the other,
# i.e., rely on an input file (and potentially include post-processing of the
# results in python)

# Heat conduction
setup_test_serial(Heat-1 MethodHeat/Heat-1  input/Heat/method_heat-1.in)
setup_test_parallel(Heat-8 MethodHeat/Heat-8  input/Heat/method_heat-8.in)

# Gravity
setup_test_serial(GravityCg-1 MethodGravity/GravityCg-1  input/Gravity/method_gravity_cg-1.in)
setup_test_parallel(GravityCg-8 MethodGravity/GravityCg-8  input/Gravity/method_gravity_cg-8.in)

# Cosmology
# Needs to be updated, see https://github.com/enzo-project/enzo-e/issues/156

# Flux correction
setup_test_serial(FluxCorrect-SMR-PPM MethodFluxCorrect/Inclined-Contact-SMR-Ppm input/FluxCorrect/inclined_contact_smr_ppm-${PREC_STRING}.in)
setup_test_serial(FluxCorrect-SMR-VL MethodFluxCorrect/Inclined-Contact-VL input/FluxCorrect/inclined_contact_smr_vl-${PREC_STRING}.in)

# PPM
setup_test_serial(PPM-1 MethodPPM/Ppm-1  input/PPM/method_ppm-1.in)
setup_test_parallel(PPM-8 MethodPPM/Ppm-8  input/PPM/method_ppm-8.in)

# VLCT
setup_test_serial_python(vlct_dual_energy_cloud vlct "input/vlct/run_dual_energy_cloud_test.py")
setup_test_serial_python(vlct_MHD_shock_tube vlct "input/vlct/run_MHD_shock_tube_test.py")
setup_test_serial_python(vlct_MHD_linear_wave vlct "input/vlct/run_MHD_linear_wave_test.py")
setup_test_serial_python(vlct_HD_linear_wave vlct "input/vlct/run_HD_linear_wave_test.py")
setup_test_serial_python(vlct_passive_advect_sound vlct "input/vlct/run_passive_advect_sound_test.py")
setup_test_parallel_python(vlct_dual_energy_shock_tube vlct "input/vlct/run_dual_energy_shock_tube_test.py")

# merge_stars
setup_test_serial_python(merge_stars_serial merge_stars/serial "input/merge_stars/run_merge_stars_test.py" "--prec=${PREC_STRING}")
setup_test_parallel_python(merge_stars_parallel merge_stars/parallel "input/merge_stars/run_merge_stars_test.py" "--prec=${PREC_STRING}")

# Convert markdown file to html file for more flexible viewing
configure_file(${CMAKE_SOURCE_DIR}/test/TEST_RESULTS.md ${PROJECT_BINARY_DIR}/TEST_RESULTS.md COPYONLY)
add_custom_target(process_test_results
                  pandoc -s -o TEST_RESULTS.html --metadata title="Enzo-E Test Results" TEST_RESULTS.md
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

