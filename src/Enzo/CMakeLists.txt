#=======================
# Define Charm++ modules
#=======================
# This section declares Charm++ module targets. These are header-only targets
# that represent the .decl.h and .def.h headers generated from the .ci
# interface files by the charm++ libraries.
#
# We also define dependencies between modules (represented by occurences of
# extern in the .ci file)
include(charm)
addCharmModule( "enzo" )
target_link_libraries(enzoCharmModule INTERFACE meshCharmModule)

# note: main_enzoCharmModule was created in a previous CMakeLists.txt file
target_link_libraries(main_enzoCharmModule INTERFACE enzoCharmModule)

#==============
# General setup
#==============

# the variable holds the commonly used root-include directory that lets files
# include headers by calling:
#    #include "Enzo/header/path.hpp"
#    #include "Cello/header/path.hpp"
# A target, called <MY-TARGET> uses of this variable in a call like:
#    target_include_directories(<MY-TARGET> PUBLIC ${ROOT_INCLUDE_DIR})
set(ROOT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../)

set(Cello_LIBS "cello_component;charm_component;control;disk;error;data;io;memory;mesh;monitor;parameters;parallel;performance;problem;compute;simulation;test_Unit")

set(External_LIBS "")
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

# Step 1: forward-declare the enzo library. This holds most of the files that
#         are turned into the enzo binary. (We will return to this shortly)
add_library(enzo STATIC charm_enzo.cpp)

# Step 2: define main_enzo
# -> this includes a machinery for driving the simulation (it's directly
#    analogous to the "tester" targets in the Cello directory)
# -> we need to define it now (and as a separate target from the binary) in
#    order to support tests that are defined in the subdirectory
add_library(main_enzo OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/../Cello/main_enzo.cpp)
# TODO: we may want to reconsider the way that the code is structured so that
# we aren't declaring targets with source files from the Cello layer

target_include_directories (main_enzo
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../Cello # for main.hpp
)
target_link_libraries(main_enzo
  PRIVATE enzo parallel cello_component monitor StdFilesystem::StdFilesystem
  PUBLIC main_enzoCharmModule # this is public because main_enzo.decl.h is in
                              # the public header file: main.hpp)
)
target_link_options(main_enzo PRIVATE ${Cello_TARGET_LINK_OPTIONS})

# Step 3: add subdirectories (turning our attention back to the enzo target)
# -> these add other libraries that will be used when linking the enzo target
# -> these also add source files to the enzo target

# add files in enzo-core to enzo target (and define a unit test)
add_subdirectory(enzo-core)

# in the following lines all introduce targets (aka libraries) that are distinct
# from the primary enzo target (later on, they are all linked to the enzo target)

# define the Enzo::assorted target
add_subdirectory(assorted)

# define the Enzo::chemistry target
add_subdirectory(chemistry)

# define the Enzo::cosmology target
add_subdirectory(cosmology)

# define the Enzo::fluid-props target
add_subdirectory(fluid-props)

# define the Enzo::gravity target
add_subdirectory(gravity)

# define the Enzo::hydro-mhd target
add_subdirectory(hydro-mhd)

# define the Enzo::initial target
add_subdirectory(initial)

# define the Enzo::inference target
add_subdirectory(inference)

# define the Enzo::io target
add_subdirectory(io)

# define the Enzo::mesh target
add_subdirectory(mesh)

# define the Enzo::particle target
add_subdirectory(particle)

# define the Enzo::tests target
# (note that the Enzo::tests target introduces problem initializers that can
# be used to setup test-problems. The BUILD_TESTING option has no impact on
# what is built from this directory. The option only affects whether binaries,
# expressly used for unit-testing are built)
add_subdirectory(tests)

# define the Enzo::utils target
add_subdirectory(utils)

# Step 4: specify additional information about the Enzo library
# (dependencies are only PUBLIC because of the way that all headers are
#  transitively included in all targets. If they weren't public, the include
#  directories would not be properly inherited)
target_link_libraries(enzo
  PUBLIC Enzo::assorted Enzo::chemistry Enzo::cosmology Enzo::fluid-props
  PUBLIC Enzo::gravity Enzo::hydro-mhd Enzo::initial Enzo::io Enzo::mesh
  PUBLIC Enzo::particle Enzo::tests Enzo::utils
  PUBLIC Enzo::inference

  PUBLIC ${Cello_LIBS} ${External_LIBS}
  PUBLIC enzoCharmModule
)
target_include_directories(enzo
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} # for enzo.hpp when linking to a target in
                                     # separate directory
)
target_link_options(enzo PRIVATE ${Cello_TARGET_LINK_OPTIONS})


# At the time of writing, cello.hpp is included in basically every source file
# (it's either directly included, or included via enzo.hpp). Thus it's okay to
# use these precompiled headers in all of the files.
target_link_libraries(enzo PRIVATE cello_component_PCH)

# At the time of wriing, enzo.hpp is included in just about every source file
# compiled as part of the enzo-target. For that reason, using
#     > if(USE_PRECOMPILED_HEADERS)
#     >    target_precompile_headers(enzo PRIVATE enzo.hpp)
#      > endif()
# actually saves a lot of extra time.
# - However, I'm VERY worried that future source files could get introduced
#   that don't include enzo.hpp. In that case, the precompiled headers would
#   allow the source file to access just about ANY in symbol in the entire
#   codebase without an appropriate include-statement.
# - While the same could still happen with the headers in cello_component_PCH,
#   I'm much less concerned since those are standard library headers and
#   charm++ headers
# - Maybe we can revisit this if we can come up with a better automatic test
#   for catching these cases (more helpful than the existing CI test). It could
#   when/if we start subdividing the enzo-target into subtargets...

#======================
# Define Binary Targets
#======================

# define the main executable
add_executable(enzo-e enzo-e.cpp)
target_link_libraries(enzo-e PRIVATE enzo main_enzo ${External_LIBS})
target_link_options(enzo-e PRIVATE ${Cello_TARGET_LINK_OPTIONS})

if (BUILD_TESTING)
  # define a unit test:
  #
  # consider removing the enzo-specific stuff from this test so that we can
  # define it entirely in the Cello layer
  add_executable(
    test_class_size "${CMAKE_CURRENT_SOURCE_DIR}/../Cello/test_class_size.cpp"
    )
  target_link_libraries(test_class_size PRIVATE enzo mesh tester_mesh)
  target_link_options(test_class_size PRIVATE ${Cello_TARGET_LINK_OPTIONS})
endif()
