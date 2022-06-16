include(charm)
addCharmModule( "mesh" )
addCharmModule( "main" )
addCharmModule( "simulation" )
addCharmModule( "main_enzo" )
addCharmModule( "main_mesh" )
addCharmModule( "main_simulation" )


function(addCelloLib LIBNAME EXTRASRCS)
  # Arguments:
  #   LIBNAME:   Name of the Cello compontent to be build.
  #   EXTRASRCS: Extra sources to be included to this target.
  # Add static library based on glob for the library name and including (potentially) excessive dependencies.

  file(GLOB SRC_FILES ${LIBNAME}*.cpp *{LIBNAME}*.hpp)
  add_library(${LIBNAME} STATIC ${SRC_FILES} ${EXTRASRCS})
  add_dependencies(${LIBNAME} simulationCharmModule mainCharmModule meshCharmModule main_meshCharmModule main_simulationCharmModule)
  target_include_directories (${LIBNAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
  target_link_options(${LIBNAME} PRIVATE ${Cello_TARGET_LINK_OPTIONS})

endfunction(addCelloLib)

addCelloLib(array "main.cpp")
target_link_libraries(array simulation)
addCelloLib(cello   "")
target_link_libraries(cello HDF5_C)
addCelloLib(charm   "")
target_link_libraries(charm HDF5_C)
addCelloLib(compute "")
target_link_libraries(compute HDF5_C)
addCelloLib(control "")
target_link_libraries(control HDF5_C)
addCelloLib(disk "main.cpp")
target_link_libraries(disk error HDF5_C cello)
addCelloLib(error "main.cpp")
target_link_libraries(error monitor Boost::filesystem)
addCelloLib(data  "main_simulation.cpp;field_face_store.F")
target_link_libraries(data HDF5_C mesh io parallel cello)
addCelloLib(io    "main_simulation.cpp")
target_link_libraries(io mesh HDF5_C)
addCelloLib(memory "main.cpp" )
target_link_libraries(memory error Boost::filesystem)
addCelloLib(mesh "main_mesh.cpp" ) 
target_link_libraries(mesh data io disk parallel cello pngwriter HDF5_C Boost::filesystem)
addCelloLib(monitor "main.cpp")
target_link_libraries(monitor parallel Boost::filesystem)
addCelloLib(parameters "main.cpp;parse.tab.c;lex.yy.c")
target_link_libraries(parameters error Boost::filesystem cello)
addCelloLib(parallel "parallel.def;main.cpp")
target_link_libraries(parallel error)
addCelloLib(performance "main.cpp")
target_link_libraries(performance memory parameters error Boost::filesystem)
if (use_papi)
  target_link_libraries(performance PAPI::papi)
endif()
addCelloLib(problem "")
target_link_libraries(problem HDF5_C)
addCelloLib(simulation "")
target_link_libraries(simulation problem HDF5_C)
addCelloLib(test_Unit "")
target_link_libraries(test_Unit parallel cello data simulation)

# Unit tests

function(addUnitTestBinary BINNAME SRCS LIBS)
  # Arguments:
  #   BINNAME:   Name of the test binary to build
  #   SRCS:      Source files
  #   LIBS:      Extra libraries required

  add_executable(${BINNAME} ${SRCS})
  target_link_libraries(${BINNAME} PUBLIC ${LIBS} test_Unit parallel cello data simulation)
  target_include_directories (${BINNAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
  target_link_options(${BINNAME} PRIVATE ${Cello_TARGET_LINK_OPTIONS})

endfunction(addUnitTestBinary)

# TODO(need help) The following commented test binaries should be included once
# someone figures out how to properly link them.
# Once fixed, the associated units tests in test/CMakeLists.txt should be uncommented.
#addUnitTestBinary(test_type "test_Type.cpp" "mesh")
#addUnitTestBinary(test_field_data "test_Field.cpp" "")
#addUnitTestBinary(test_field_descr "test_FieldDescr.cpp" "")
#addUnitTestBinary(test_field_face "test_FieldFace.cpp" "")
#addUnitTestBinary(test_grouping "test_Grouping.cpp" "")
#addUnitTestBinary(test_itindex "test_ItIndex.cpp" "")
#addUnitTestBinary(test_scalar "test_Scalar.cpp" "")
#addUnitTestBinary(test_enzo_units "test_EnzoUnits.cpp" "")
addUnitTestBinary(test_carr_collec "test_CArrCollec.cpp" "array")
addUnitTestBinary(test_cello_array "test_CelloArray.cpp" "array")
addUnitTestBinary(
  test_string_ind_rd_only_map "test_StringIndRdOnlyMap.cpp" "array"
)
#addUnitTestBinary(test_colormap "test_Colormap.cpp" "io")
addUnitTestBinary(test_error "test_Error.cpp" "error")
#addUnitTestBinary(test_field "test_Field.cpp" "")
addUnitTestBinary(test_memory "test_Memory.cpp" "memory")
addUnitTestBinary(test_monitor "test_Monitor.cpp" "monitor")
#addUnitTestBinary(test_particle "test_Particle.cpp" "")
#addUnitTestBinary(test_ "test_.cpp" "")
