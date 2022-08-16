
function(setup_test_dir TESTDIR)
  # make test directory within build tree (will be used for test output/current working directory)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  # link in input directory because of fixed input file hierarchy 
  file(CREATE_LINK ${PROJECT_SOURCE_DIR}/input ${PROJECT_BINARY_DIR}/test/${TESTDIR}/input SYMBOLIC)
endfunction()

# Functions for setting up C++-based tests
#-----------------------------------------
# these tests rely upon the testing machinery defined in Cello.
# - These tests effectively just execute C++ binaries (they're either dedicated
#   unit tests or the Enzo-E binary)
# - At the moment, the exit codes of these tests don't always report the
#   correct exit code. Therefore, a lightweight python script is used to
#   actually execute the binary, parse the results from stdout, and then
#   provide the correct exit code. This python script has the following path:
set(CPP_TEST_RUNNER ${PROJECT_SOURCE_DIR}/tools/run_cpp_test.py)

# Function that sets up a unit test (each unit test is organized into a
# separate C++ binary)
function(setup_test_unit TESTNAME TESTDIR TESTBIN)
  setup_test_dir(${TESTDIR})
  set(FULLTESTDIR ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND python3 ${CPP_TEST_RUNNER} --output-dump ${FULLTESTDIR}/${TESTNAME}.log $<TARGET_FILE:${TESTBIN}>
    WORKING_DIRECTORY ${FULLTESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;unit" )
endfunction()

# Function that sets up a test that directly calls Enzo-E using a single
# compute unit. These tests pass or fail based on whether they run to
# completion and whether expectations about the stopping time or cycle
# (expressed in the parameter file) are satisfied
function(setup_test_serial TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  set(FULLTESTDIR ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND python3 ${CPP_TEST_RUNNER} --output-dump ${FULLTESTDIR}/${TESTNAME}.log $<TARGET_FILE:enzo-e> ${INFILE}
    WORKING_DIRECTORY ${FULLTESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;enzo" )
endfunction()

# Just like setup_test_serial except that Enzo-E is executed with multiple
# compute units
function(setup_test_parallel TESTNAME TESTDIR INFILE)
  setup_test_dir(${TESTDIR})
  set(FULLTESTDIR ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND python3 ${CPP_TEST_RUNNER} --output-dump ${FULLTESTDIR}/${TESTNAME}.log ${PARALLEL_LAUNCHER} ${PARALLEL_LAUNCHER_NPROC_ARG} ${PARALLEL_LAUNCHER_NPROC} $<TARGET_FILE:enzo-e> ${INFILE}
    WORKING_DIRECTORY ${FULLTESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel;enzo" )
endfunction()

# Functions setting up YT_BASED_TESTS
# -----------------------------------
# - these tests execute a python script that both executes the enzo-e binary
#   and then performs analysis on the outputs

# these tests all require that python has yt, libconf, and h5py installed.
# - We express this requirement in a "fixture" called YT_DEPEND
# - all tests that require the YT_DEPEND fixture will only be run if the
#   has_yt_dependencies test passes. Otherwise, they get skipped.

# a test for whether the python interpreter has necessary dependencies
add_test(NAME python_has_yt_and_other_modules
  COMMAND python3 -c "import yt; import libconf; import h5py")
set_tests_properties(python_has_yt_and_other_modules
  PROPERTIES FIXTURES_SETUP YT_DEPEND)

# Serial (single compute element) test wrapper functions calling a script (e.g.
# to include output analysis). All optional argment (i.e. the ones not being
# named) are automatically passed to the command being called.
function(setup_test_serial_python TESTNAME TESTDIR TESTSCRIPT)
  setup_test_dir(${TESTDIR})
  add_test(
    NAME ${TESTNAME}
    COMMAND python3 ${TESTSCRIPT} --launch_cmd=${PROJECT_BINARY_DIR}/bin/enzo-e ${ARGN}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "serial;enzo;regression"
    FIXTURES_REQUIRED YT_DEPEND)
endfunction()

# Parallel (using multiple compute elements) test wrapper functions calling a
# script (e.g., to include output analysis)
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
    COMMAND python3 ${TESTSCRIPT} --launch_cmd ${LAUNCH_CMD} ${ARGN}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TESTDIR})
  set_tests_properties(${TESTNAME} PROPERTIES LABELS "parallel;enzo;regression"
    FIXTURES_REQUIRED YT_DEPEND)
endfunction()

############################### UNIT TESTS ####################################
# The following tests are self-contained binary of individual components

setup_test_unit(CArrCollec ArrayComponent/CArrCollec test_carr_collec)
setup_test_unit(Array ArrayComponent/Array test_cello_array)
setup_test_unit(StringIndRdOnlyMap ArrayComponent/StringIndRdOnlyMap test_string_ind_rd_only_map)

setup_test_unit(CelloType Cello/Type test_type)

setup_test_unit(FileHDF5 DiskComponent/FileHDF5 test_file_hdf5)

setup_test_unit(Error ErrorComponent/Error test_error)

setup_test_unit(Colormap IOComponent/Colormap test_colormap)
setup_test_unit(Schedule IOComponent/Schedule test_schedule)

setup_test_unit(Data-Particle DataComponent/Particle test_particle)
setup_test_unit(Data-Scalar DataComponent/Scalar test_scalar)
setup_test_unit(Data-Field-Data DataComponent/FieldData test_field_data)
setup_test_unit(Data-Field-Descr DataComponent/FieldDescr test_field_descr)
setup_test_unit(Data-Field DataComponent/Field test_field)
setup_test_unit(Data-Field-Face DataComponent/FieldFace test_field_face)
setup_test_unit(Data-Grouping DataComponent/Grouping test_grouping)
setup_test_unit(Data-ItIndex DataComponent/ItIndex test_itindex)

setup_test_unit(Memory MemoryComponent/Memory test_memory)

setup_test_unit(Monitor MonitorComponent/Monitor test_monitor)

setup_test_unit(
  Parameters-Parameters ParametersComponent/Parameters test_parameters
)
# we need to pass an argument to the following test (a path to a config file)
#setup_test_unit(Parameters-Parse ParametersComponent/Parse test_parse)

setup_test_unit(
  Performance-Performance PerformanceComponent/Performance test_performance
)
setup_test_unit(Performance-Timer PerformanceComponent/Timer test_timer)
if (use_papi)
  setup_test_unit(Performance-Papi PerformanceComponent/Papi test_papi)
endif()

setup_test_unit(EnzoUnits UnitsComponent/EnzoUnits test_enzo_units)

# TODO: sort the following test by component
setup_test_unit(Assorted-class_size Assorted/class_size test_class_size)
setup_test_unit(Assorted-Data Assorted/Data test_data)
setup_test_unit(Assorted-ItFace Assorted/ItFace test_it_face)
setup_test_unit(Assorted-ItChild Assorted/ItChild test_it_child)
setup_test_unit(Assorted-BlockTrace Assorted/BlockTrace test_block_trace)
setup_test_unit(Assorted-Sync Assorted/Sync test_sync)
setup_test_unit(Assorted-Face Assorted/Face test_face)
setup_test_unit(Assorted-FaceFluxes Assorted/FaceFluxes test_face_fluxes)
setup_test_unit(Assorted-FluxData Assorted/FluxData test_flux_data)
setup_test_unit(
  Assorted-ProlongLinear Assorted/ProlongLinear test_prolong_linear
)
setup_test_unit(Assorted-Refresh Assorted/Refresh test_refresh)
setup_test_unit(Assorted-Mask Assorted/Mask test_mask)
setup_test_unit(Assorted-Value Assorted/Value test_value)
setup_test_unit(Assorted-Box Assorted/Box test_box)

# TODO we need to fix the following tests (see commented units tests in
# src/Cello/CMakeLists.txt)
#
# the test runs, but one of the checks fails.
#setup_test_unit(Assorted-Adapt Assorted/Adapt test_adapt)
#
# the following test is broken (an API change is not reflected in the tests)
#setup_test_unit(Assorted-Index Assorted/Index test_index)


#setup_test_unit( Component/ test_)

############################### ENZO-E TESTS ##################################
# The following tests will call the enzo-e binary in one way or the other,
# i.e., rely on an input file (and potentially include post-processing of the
# results in python)

# AMR PPM Adapt
setup_test_parallel(AmrPpm AmrPpm/Adapt  input/Adapt/adapt-L5-P1.in)

# Boundaries
setup_test_parallel(Bound-Reflect-2D BoundaryConditions/Reflecting-2D  input/Boundary/boundary_reflecting-2d.in)
setup_test_parallel(Bound-Periodic-2D BoundaryConditions/Periodic-2D  input/Boundary/boundary_periodic-2d.in)
setup_test_parallel(Bound-Outflow-2D BoundaryConditions/Outflow-2D  input/Boundary/boundary_outflow-2d.in)
setup_test_parallel(Bound-Reflect-3D BoundaryConditions/Reflecting-3D  input/Boundary/boundary_reflecting-3d.in)
setup_test_parallel(Bound-Periodic-3D BoundaryConditions/Periodic-3D  input/Boundary/boundary_periodic-3d.in)
setup_test_parallel(Bound-Outflow-3D BoundaryConditions/Outflow-3D  input/Boundary/boundary_outflow-3d.in)

# Cosmology
# Needs to be updated, see https://github.com/enzo-project/enzo-e/issues/156

# Flux correction
setup_test_serial(FluxCorrect-SMR-PPM MethodFluxCorrect/Inclined-Contact-SMR-Ppm input/FluxCorrect/inclined_contact_smr_ppm-${PREC_STRING}.in)
setup_test_serial(FluxCorrect-SMR-VL MethodFluxCorrect/Inclined-Contact-VL input/FluxCorrect/inclined_contact_smr_vl-${PREC_STRING}.in)

# Isolated galaxy
setup_test_serial(GasDisk IsolatedGalaxy/GasDisk  input/IsolatedGalaxy/method_isolatedgalaxy.in)
#setup_test_serial(GasDisk-Halo IsolatedGalaxy/GasDisk-Halo  input/IsolatedGalaxy/method_isolatedgalaxy-particles.in)

# Gravity
setup_test_serial(GravityCg-1 MethodGravity/GravityCg-1  input/Gravity/method_gravity_cg-1.in)
setup_test_parallel(GravityCg-8 MethodGravity/GravityCg-8  input/Gravity/method_gravity_cg-8.in)

# Heat conduction
setup_test_serial(Heat-1 MethodHeat/Heat-1  input/Heat/method_heat-1.in)
setup_test_parallel(Heat-8 MethodHeat/Heat-8  input/Heat/method_heat-8.in)

# Initial
setup_test_serial(Music-111 InitialMusic/Music-111  input/InitialMusic/initial_music-111.in)
setup_test_serial(Music-222 InitialMusic/Music-222  input/InitialMusic/initial_music-222.in)
setup_test_serial(Music-444 InitialMusic/Music-444  input/InitialMusic/initial_music-444.in)
setup_test_serial(Music-211 InitialMusic/Music-211  input/InitialMusic/initial_music-211.in)
setup_test_serial(Music-121 InitialMusic/Music-121  input/InitialMusic/initial_music-121.in)
setup_test_serial(Music-112 InitialMusic/Music-112  input/InitialMusic/initial_music-112.in)
setup_test_serial(Music-411 InitialMusic/Music-411  input/InitialMusic/initial_music-411.in)
setup_test_serial(Music-141 InitialMusic/Music-141  input/InitialMusic/initial_music-141.in)
setup_test_serial(Music-114 InitialMusic/Music-114  input/InitialMusic/initial_music-114.in)

# Output
setup_test_parallel(Output-Stride-1 Output/Output-Stride-1  input/Output/output-stride-1.in)
setup_test_parallel(Output-Stride-2 Output/Output-Stride-2  input/Output/output-stride-2.in)
setup_test_parallel(Output-Stride-4 Output/Output-Stride-4  input/Output/output-stride-4.in)
setup_test_parallel(Output-Headers  Output/Output-Headers   input/Output/output-headers.in)

# Particles
setup_test_parallel(Particle-X  Particles/X   input/Particle/test_particle-x.in)
setup_test_parallel(Particle-Y  Particles/Y   input/Particle/test_particle-y.in)
setup_test_parallel(Particle-XY  Particles/XY   input/Particle/test_particle-xy.in)
setup_test_parallel(Particle-Circle  Particles/Circle   input/Particle/test_particle-circle.in)
setup_test_parallel(Particle-AMR-Static Particles/AMR-Static   input/Particle/test_particle-amr-static.in)
setup_test_parallel(Particle-AMR-Dynamic Particles/AMR-Dynamic   input/Particle/test_particle-amr-dynamic.in)

# Performance
setup_test_parallel(Performance-Initial-PNG  Performance/InitialPng   input/HelloWorld/initial_png.in)

# PPM
setup_test_serial(PPM-1 MethodPPM/Ppm-1  input/PPM/method_ppm-1.in)
setup_test_parallel(PPM-8 MethodPPM/Ppm-8  input/PPM/method_ppm-8.in)
setup_test_serial(PPM-1_color MethodPPM/Ppm-1_color  input/PPM/method_ppm-1_color.in)


# define yt-based tests

# VLCT
setup_test_serial_python(vlct_dual_energy_cloud vlct "input/vlct/run_dual_energy_cloud_test.py")
setup_test_serial_python(vlct_MHD_shock_tube vlct "input/vlct/run_MHD_shock_tube_test.py")
setup_test_serial_python(vlct_MHD_linear_wave vlct "input/vlct/run_MHD_linear_wave_test.py")
setup_test_serial_python(vlct_HD_linear_wave vlct "input/vlct/run_HD_linear_wave_test.py")
setup_test_serial_python(vlct_passive_advect_sound vlct "input/vlct/run_passive_advect_sound_test.py")
setup_test_parallel_python(vlct_dual_energy_shock_tube vlct "input/vlct/run_dual_energy_shock_tube_test.py")

# Gravity (with VLCT)
setup_test_serial_python(gravity_vlct_stable_Jeans_wave gravity "input/Gravity/run_stable_jeans_wave_test.py")

# merge_sinks
setup_test_serial_python(merge_sinks_stationary_serial merge_sinks/stationary/serial "input/merge_sinks/run_merge_sinks_test.py" "--prec=${PREC_STRING}" "--ics_type=stationary")
setup_test_parallel_python(merge_sinks_stationary_parallel merge_sinks/stationary/parallel "input/merge_sinks/run_merge_sinks_test.py" "--prec=${PREC_STRING}" "--ics_type=stationary")
setup_test_serial_python(merge_sinks_drift_serial merge_sinks/drift/serial "input/merge_sinks/run_merge_sinks_test.py" "--prec=${PREC_STRING}" "--ics_type=drift")
setup_test_parallel_python(merge_sinks_drift_parallel merge_sinks/drift/parallel "input/merge_sinks/run_merge_sinks_test.py" "--prec=${PREC_STRING}" "--ics_type=drift")

# accretion
setup_test_serial_python(threshold_accretion_serial accretion/threshold/serial "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=threshold")
setup_test_parallel_python(threshold_accretion_parallel accretion/threshold/parallel "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=threshold")
setup_test_serial_python(bondi_hoyle_accretion_serial accretion/bondi_hoyle/serial "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=bondi_hoyle")
setup_test_parallel_python(bondi_hoyle_accretion_parallel accretion/bondi_hoyle/parallel "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=bondi_hoyle")
setup_test_serial_python(flux_accretion_serial accretion/flux/serial "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=flux")
setup_test_parallel_python(flux_accretion_parallel accretion/flux/parallel "input/accretion/run_accretion_test.py" "--prec=${PREC_STRING}" "--flavor=flux")

# Shu Collapse
setup_test_parallel_python(shu_collapse shu_collapse "input/shu_collapse/run_shu_collapse.py" "--prec=${PREC_STRING}")

# BB Test
setup_test_parallel_python(bb_test bb_test "input/bb_test/run_bb_test.py" "--prec=${PREC_STRING}")

# Grackle tests
#
# because Grackle is a separate library, the user needs to specify the path
# to the data directory where the grackle data files are installed. if this is
# not specified, the Grackle tests are skipped

if(NOT (DEFINED GRACKLE_INPUT_DATA_DIR))
  set(GRACKLE_INPUT_DATA_DIR "")
endif()

if(USE_GRACKLE AND (GRACKLE_INPUT_DATA_DIR STREQUAL ""))
  message(STATUS
    "   IMPORTANT: No tests involving Grackle have been defined (even though Grackle is being used). To enable these tests, try `-DGRACKLE_INPUT_DATA_DIR=/path/to/grackle/data/dir`."
  )
elseif(USE_GRACKLE)
  # (non-yt-based test)
  setup_test_parallel_python(grackle_cooling_dt grackle
    "input/Grackle/run_endtime_grackle_test.py"
    "--grackle-data-file=${GRACKLE_INPUT_DATA_DIR}/CloudyData_UVB=HM2012_shielded.h5"
    "--nominal-config-path=input/Grackle/method_grackle_cooling_dt.in"
    # by using relative paths in the following 2 arguments, we effectively
    # specify output files in the testing directory
    "--generate-config-path=./temp_grackle_cooling_dt.in"
    "--output-dump=./grackle_cooling_dt.log"
  )

  setup_test_parallel_python(grackle_general grackle
    "input/Grackle/run_general_grackle_test.py"
    "--prec=${PREC_STRING}"
    "--grackle-data-file=${GRACKLE_INPUT_DATA_DIR}/CloudyData_UVB=HM2012_shielded.h5"
    # by using relative paths in the following 2 arguments, we effectively
    # specify output files in the testing directory
    "--generate-config-path=./temp_grackle_general.in"
    "--output-dump=./grackle_general_test.log"
  )
endif()

# Convert markdown file to html file for more flexible viewing
configure_file(${CMAKE_SOURCE_DIR}/test/TEST_RESULTS.md ${PROJECT_BINARY_DIR}/TEST_RESULTS.md COPYONLY)
add_custom_target(process_test_results
                  pandoc -s -o TEST_RESULTS.html --metadata title="Enzo-E Test Results" TEST_RESULTS.md
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

