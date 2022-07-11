
include(charm)
addCharmModule( "enzo" )

file(GLOB SRC_FILES enzo*.cpp *F *enzo*hpp fortran.h fortran_types.h FofLib.cpp)

# TODO(pgrete) Check if that file can be removed for good. It also not being used in the SCons build.
# -> PPML_sub.F:(.text+0x5ebf): multiple definition of `qdd6_'; CMakeFiles/enzo-e.dir/PPML_QDD6.F.o:PPML_QDD6.F:(.text+0x0): first defined here
list(FILTER SRC_FILES EXCLUDE REGEX "PPML_QDD6")

set(Cello_LIBS "cello_component;charm_component;control;disk;error;data;io;memory;mesh;monitor;parameters;parallel;performance;problem;compute;simulation;test_Unit")

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

#=======================
# Define Library Targets
#=======================
add_subdirectory(EnzoRiemann) # define the riemann library


# build an enzo library that is then turned into an enzo-e binary
add_library(enzo STATIC ${SRC_FILES})
add_dependencies(enzo enzoCharmModule)
# dependencies are only PUBLIC because of the way that all headers are
# transitively included in all targets. If they weren't public, the include
# directories would not be properly inherited
target_link_libraries(enzo PUBLIC riemann ${Cello_LIBS} ${External_LIBS})
target_include_directories(enzo
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR} # enzo.decl.h is in a public header file
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} # for enzo.hpp when linking to a target in
                                     # separate directory
)
target_link_options(enzo PRIVATE ${Cello_TARGET_LINK_OPTIONS})


# this includes a lot of the machinery for driving the simulation (it's
# directly analogous to the "tester" targets in the Cello directory)
#
# we define this as a separate target from the binary in order to support tests
add_library(main_enzo OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/../Cello/main_enzo.cpp)
add_dependencies(main_enzo main_enzoCharmModule)
# TODO: we may want to reconsider the way that the code is structured so that
# we aren't declaring targets with source files from the other layer
target_include_directories (main_enzo
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../Cello # for main_enzo.decl.h, (which is
                                              # in a public header: main.hpp)
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../Cello # for main.hpp
)
target_link_libraries(main_enzo
  PRIVATE enzo parallel cello_component monitor Boost::filesystem
)
target_link_options(main_enzo PRIVATE ${Cello_TARGET_LINK_OPTIONS})

#======================
# Define Binary Targets
#======================

# define the main executable
add_executable(enzo-e enzo-e.cpp)
target_link_libraries(enzo-e PRIVATE enzo main_enzo ${External_LIBS})
target_link_options(enzo-e PRIVATE ${Cello_TARGET_LINK_OPTIONS})

# define some unit tests:

# This tests code
add_executable(test_enzo_units "test_EnzoUnits.cpp")
target_link_libraries(test_enzo_units PRIVATE enzo main_enzo)
target_link_options(test_enzo_units PRIVATE ${Cello_TARGET_LINK_OPTIONS})

# consider removing the enzo-specific stuff from this test so that we can
# define it entirely in the Cello layer
add_executable(
  test_class_size "${CMAKE_CURRENT_SOURCE_DIR}/../Cello/test_class_size.cpp"
)
target_link_libraries(test_class_size PRIVATE enzo mesh tester_mesh)
target_link_options(test_class_size PRIVATE ${Cello_TARGET_LINK_OPTIONS})
