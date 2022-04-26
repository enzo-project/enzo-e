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

addCelloLib(cello   "")
target_link_libraries(cello HDF5_C)
addCelloLib(charm   "")
target_link_libraries(charm HDF5_C)
addCelloLib(control "")
target_link_libraries(control HDF5_C)
addCelloLib(disk "main.cpp")
target_link_libraries(disk HDF5_C)
addCelloLib(error "main.cpp")
addCelloLib(data  "main_simulation.cpp;field_face_store.F")
target_link_libraries(data HDF5_C)
addCelloLib(io    "main_simulation.cpp")
target_link_libraries(io HDF5_C)
# need to figure out what to do with enzo from main_enzo in this src tree
#addCelloLib(main  "")
addCelloLib(memory "main.cpp" ) 
addCelloLib(mesh "main_mesh.cpp" ) 
target_link_libraries(mesh pngwriter HDF5_C)
addCelloLib(monitor "main.cpp")
addCelloLib(parameters "main.cpp;parse.tab.c;lex.yy.c")
addCelloLib(parallel "parallel.def")
addCelloLib(performance "main.cpp")
if (use_papi)
  target_link_libraries(performance PAPI::papi)
endif()
addCelloLib(problem "")
target_link_libraries(problem HDF5_C)
addCelloLib(compute "")
target_link_libraries(compute HDF5_C)
addCelloLib(simulation "")
target_link_libraries(simulation HDF5_C)
addCelloLib(test_Unit "")

# Unit tests

function(addUnitTestBinary BINNAME SRCS LIBS)
  # Arguments:
  #   BINNAME:   Name of the test binary to build
  #   SRCS:      Source files
  #   LIBS: Extra libraries required

  add_executable(${BINNAME} ${SRCS})
  target_link_libraries(${BINNAME} PUBLIC ${LIBS} test_Unit cello parallel error monitor Boost::filesystem)
  target_include_directories (${BINNAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
  target_link_options(${BINNAME} PRIVATE ${Cello_TARGET_LINK_OPTIONS})

endfunction(addUnitTestBinary)

addUnitTestBinary(test_memory "test_Memory.cpp" "memory")
addUnitTestBinary(test_monitor "test_Monitor.cpp" "monitor")
#addUnitTestBinary(test_cello_array "test_CelloArray.cpp;array.hpp;array_CelloArray.hpp" "mesh")
