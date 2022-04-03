
include(charm)
addCharmModule( "enzo" )

file(GLOB SRC_FILES enzo*.cpp *F *enzo*hpp fortran.h fortran_types.h FofLib.cpp)

# TODO(pgrete) Check if that file can be removed for good. It also not being used in the SCons build.
# -> PPML_sub.F:(.text+0x5ebf): multiple definition of `qdd6_'; CMakeFiles/enzo-e.dir/PPML_QDD6.F.o:PPML_QDD6.F:(.text+0x0): first defined here
list(FILTER SRC_FILES EXCLUDE REGEX "PPML_QDD6")

set(Cello_LIBS "cello;charm;control;disk;error;data;io;memory;mesh;monitor;parameters;parallel;performance;problem;compute;simulation;test_Unit")

set(External_LIBS "pngwriter;Boost::filesystem;PNG::PNG;HDF5_C")
if (USE_GRACKLE)
  LIST(APPEND External_LIBS "Grackle::Grackle")
endif()
if (use_papi)
  LIST(APPEND External_LIBS "PAPI::papi")
endif()
if (use_jemalloc)
  LIST(APPEND External_LIBS "jemalloc::jemalloc")
endif()

add_subdirectory(EnzoRiemann) # define the riemann library

# Following the SCons build and first build an enzo library that is then turned into an enzo-e binary
# (Note, the first attempt to directly build the binary resulted in failed vtable lookups in GrackleMethod and GrackleInputTest classes)
add_library(enzo STATIC ${SRC_FILES})
add_dependencies(enzo enzoCharmModule main_enzoCharmModule simulationCharmModule)

target_link_libraries(enzo
  PUBLIC riemann # this is public so that dependents of the enzo library will
                 # inherit the necessary include directories (only necessary
                 # b/c the riemann library defines a charm++ module)
  PRIVATE ${Cello_LIBS} ${External_LIBS}
)
target_include_directories (enzo PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
target_link_options(enzo PRIVATE ${Cello_TARGET_LINK_OPTIONS})

add_executable(enzo-e enzo-e.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../Cello/main_enzo.cpp)
add_dependencies(enzo-e enzoCharmModule main_enzoCharmModule simulationCharmModule)
target_link_libraries(enzo-e PRIVATE enzo ${External_LIBS})
# TOOD(pgrete) need a cleaner place for the generated headers (current in Cello folder) or export header location
target_include_directories (enzo-e PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CHARM_INCLUDE_DIRS})
target_link_options(enzo-e PRIVATE ${Cello_TARGET_LINK_OPTIONS})

