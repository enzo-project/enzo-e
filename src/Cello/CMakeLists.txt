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
addCelloLib(charm   "")
addCelloLib(control "")
addCelloLib(disk "main.cpp")
addCelloLib(error "main.cpp")
addCelloLib(data  "main_simulation.cpp;field_face_store.F")
addCelloLib(io    "main_simulation.cpp")
# need to figure out what to do with enzo from main_enzo in this src tree
#addCelloLib(main  "")
addCelloLib(memory "main.cpp" ) 
addCelloLib(mesh "main_mesh.cpp" ) 
target_link_libraries(mesh pngwriter)
addCelloLib(monitor "main.cpp")
addCelloLib(parameters "main.cpp;parse.tab.c;lex.yy.c")
addCelloLib(parallel "parallel.def")
addCelloLib(performance "main.cpp")
if (use_papi)
  target_link_libraries(performance PAPI::papi)
endif()
addCelloLib(problem "")
addCelloLib(compute "")
addCelloLib(simulation "")
addCelloLib(test_Unit "")


add_executable(test_monitor test_Monitor.cpp)
target_link_libraries(test_monitor PUBLIC monitor test_Unit cello parallel error Boost::filesystem)
target_include_directories (test_monitor PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
target_link_options(test_monitor PRIVATE ${Cello_TARGET_LINK_OPTIONS})

