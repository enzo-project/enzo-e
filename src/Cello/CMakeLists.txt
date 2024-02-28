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
addCharmModule( "mesh" )
addCharmModule( "main" )

addCharmModule( "simulation" )
target_link_libraries(simulationCharmModule INTERFACE meshCharmModule)

addCharmModule( "main_mesh" )
target_link_libraries(main_meshCharmModule INTERFACE meshCharmModule)

addCharmModule( "main_simulation" )
target_link_libraries(main_simulationCharmModule
  INTERFACE meshCharmModule simulationCharmModule
)

addCharmModule( "main_enzo" )
target_link_libraries(main_enzoCharmModule
  INTERFACE meshCharmModule simulationCharmModule
) # we append dependency on enzoCharmModule later (when it's defined)

#==============
# General setup
#==============

# In the future, it would be worth refactoring things so that every target
# doesn't depend on every charm module. For example, the error target probably
# doesn't need to depend on everything else,

function(setupCelloLibHelper LIBNAME)
  # Arguments:
  #   LIBNAME: Name of the Cello library that metadata should be modified for
  # This is just a helper function that was split off of addCelloLib for use
  # with a handful of special targets.

  target_link_options(${LIBNAME} PRIVATE ${Cello_TARGET_LINK_OPTIONS})

  # the following command both:
  # - introduces dependencies onto .ci files used to generate .decl.h and
  #   .def.h header files
  # - it also implicitly handles include directories (so that the .decl.h and
  #   .def.h header files can be found)
  # These are listed as PUBLIC dependencies for now, for safety. But in the
  # future (when we pair down the list of dependencies), we can probably change
  # the specifier (to INTERFACE/PRIVATE) based on whether the generated headers
  # are in public headers
  target_link_libraries(${LIBNAME}
    PUBLIC simulationCharmModule mainCharmModule meshCharmModule main_meshCharmModule main_simulationCharmModule
  )

  # the following statment ensures that targets (that depend on ${LIBNAME} and
  # are defined in other directories) will look for header files in the
  # in correct directory.
  target_include_directories (${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

endfunction(setupCelloLibHelper)

function(addCelloLib LIBNAME EXTRASRCS)
  # Arguments:
  #   LIBNAME:   Name of the Cello compontent to be build.
  #   EXTRASRCS: Extra sources to be included to this target.
  # Add static library based on glob for the library name and including
  # (potentially) excessive dependencies.

  # Get the list of source files
  # - we do this using GLOB patterns. This approach is not recommended by the
  #   authors of CMake (their recommendation is to explicitly list all files
  #   that must be compiled).
  # - Some of the disadvantages of this approach are mitigated by inclusion of
  #   the CONFIGURE_DEPENDS flag.
  # - See the CMake Primer section of the documentation for more details
  file(GLOB SRC_FILES ${LIBNAME}*.cpp *{LIBNAME}*.hpp)
  add_library(${LIBNAME} STATIC ${SRC_FILES} ${EXTRASRCS})
  setupCelloLibHelper(${LIBNAME})

endfunction(addCelloLib)

message("")
message("link options: ${Cello_TARGET_LINK_OPTIONS}")


# NOTE: parallel and charm_component both reference mesh_Index.hpp

#===============================================
# Defining Cello-Component Targets step 1 (of 3)
#===============================================
# This section includes placeholder target(s) that are effectively header-only
#
# Prior to CMake version 3.23, we can't actually specify the header files for
# a header-only interface-library. Because of this and the fact that all Cello
# files are currently organized in a flat directory, we don't actuall need to
# define the targets in this section. However, we make targets for these
# components anyway for explicitness and in case we ever reorganize or
# transform these libraries so that they are no-longer header-only

add_library(parallel INTERFACE)
target_link_libraries(parallel
  INTERFACE cello_component charm_component error
  INTERFACE CHARM::CHARM_HEADERS # charm++.hpp included in a public header file
)

add_library(view INTERFACE)
target_link_libraries(view INTERFACE error monitor)


# A brief segue: define a common target, cello_component_PCH, that is used to
# propagate a precompiled header
# - if USE_PRECOMPILED_HEADERS is OFF, then linking against this target does
#   nothing!
# - A precompiled header stores the contents of preparsed headers in a binary
#   format. They can significantly reduce compilation times when they are used
#   to define large headers that are included in MANY files.
# - NOTE: there's still a compile-time cost to including unnecessary headers
#   (For example, even when the <charm++.h> header is precompiled, a bunch of
#   time is still spent instantiate functions)
#
# Using a precompiled header will forcibly include all preparsed headers
# inside a source file. If care isn't taken, this can unintentionally lead to
# a source file knowing about symbols that it hasn't properly included.
# Obviously, we REALLY want to avoid that scenario.
#
# For that reason, we are going to do something simple, we are just going to
# precompile some of the expensive external headers included by cello.hpp.
# - we will only use this for targets for which ALL source files include the
#   cello.hpp file (as a rule of thumb, this should just be used when the
#   cello_component target, defined below, is a public dependency). This saves
#   time parsing the <vector> and <charm++.h> headers.
# - In the worst case scenario, we are limiting the symbols can know about
#   (that it shouldn't) to stuff in the std:: namespace and stuff from charm++
#
# Because it's good practice to make precompiled headers opt-in (since a given
# file can only use 1 precompiled header), all targets that link to this to
# target use a PRIVATE dependence (to avoid transitive dependence)
add_library(cello_component_PCH INTERFACE)
if(USE_PRECOMPILED_HEADERS)
  target_precompile_headers(cello_component_PCH INTERFACE
    <random> # this takes a shockingly long time to compile (and really
             # shouldn't be in cello.hpp)
    <string> <vector> <charm++.h> "\"pup_stl.h\""
  )
endif()

#===============================================
# Defining Cello-Component Targets step 2 (of 3)
#===============================================
# This section involves the definition of targets that actually correspond to
# concrete static libraries.
#
# To determine the link_libraries we have followed the following procedure for
# each of the targets in this section:
#   - For each target, we have gone through and declared PUBLIC dependencies
#     based on the include statements in the public header file. For example,
#     for the compute component, we list all of the components in the
#     compute.hpp header file.
#   - Because the codebase is currently organized such that all headers for a
#     given component are included in all source files in any dependent
#     component, we also had to selectively declare some other PUBLIC
#     dependencies (some light refactoring would let us clean this up -
#     especially when it comes to external dependencies)
#   - For a given target, if headers for a dependency only show up in a source
#     file of the target (and not one of the target's headers), we can then
#     declare the dependency to be PRIVATE
#
# Note: Some components seem to include the headers of some unnecessary
# dependent targets (e.g. does the compute component NEED to directly depend on
# the disk component?). We could definitely remove some of these dependencies
# (and rely upon on CMake's ability to handle such transitive dependencies for
# us). This might help speed up incremental re-compilation


# first, we define 2 targets that have "special names"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These names were chosen to minimize confusion (we have other targets that
# share similar names)


# note: in the near future, we may make a separate target called something like
# CELLO_COMBINED to represent the complete library that is composed of the
# components in this directory
file(GLOB CELLO_COMP_SRC_FILES cello*.cpp *cello*.hpp)
add_library(cello_component ${CELLO_COMP_SRC_FILES})
setupCelloLibHelper(cello_component)
target_link_libraries (cello_component
  PUBLIC CHARM::CHARM_HEADERS # charm++.h included in public header
  PRIVATE Boost::filesystem
  PRIVATE charm_component error simulation # these private dependencies arise
                                           # from source files
  PRIVATE cello_component_PCH # cello component-precompiled header
)



file(GLOB CHARM_COMP_SRC_FILES charm*.cpp *charm*.hpp)
add_library(charm_component ${CHARM_COMP_SRC_FILES})
setupCelloLibHelper(charm_component)
target_link_libraries (charm_component
  PUBLIC error
  PUBLIC CHARM::CHARM_HEADERS # charm++.h included in public header
  PRIVATE mesh simulation # these private dependencies arise from source files
)
target_include_directories (charm_component
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR} # simulation.decl.h included in public
                                     # header
)

# now, we'll define the remaining targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# these targets all correspond to files that are named in a similar way which
# makes auto-generation of these libraries easier.

addCelloLib(compute "")
target_link_libraries(compute
  PUBLIC cello_component charm_component error performance monitor memory
         parallel disk io parameters problem mesh view data simulation
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


addCelloLib(control "")
target_link_libraries(control
  PUBLIC cello_component charm_component
  PRIVATE simulation mesh # these arise from source files
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


# this probably doesn't need to depend on things like disk or compute
addCelloLib(data  "field_face_store.F")
target_link_libraries(data
  PUBLIC cello_component charm_component error memory performance monitor
         parallel disk problem mesh parameters io compute simulation view
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


addCelloLib(disk "")
target_link_libraries(disk
  PRIVATE HDF5_C
  PRIVATE pngwriter
  PUBLIC error
  )
# the disk_utils functions requires that all linked libraries are all "linked"
# against the view library (the disk library itself does not need to link
# to the view library because none of the symbols defined in that library are
# used within the disk library's source files)
target_link_libraries(disk INTERFACE view)


addCelloLib(error "")
target_link_libraries(error
  PUBLIC CHARM::CHARM_HEADERS # charm++.h included in public header file
  PUBLIC performance monitor parallel
)


addCelloLib(io "")
target_link_libraries(io
  PUBLIC cello_component charm_component parallel monitor performance error
         parameters disk problem compute mesh view data simulation
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


addCelloLib(memory "")
target_link_libraries(memory PUBLIC error parallel monitor)


addCelloLib(mesh "")
target_link_libraries(mesh
  PUBLIC cello_component charm_component performance monitor parallel disk
         parameters control io problem compute simulation view data
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


addCelloLib(monitor "")
# the next line is for including auto_config.def
target_include_directories (monitor PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(monitor
  PUBLIC CHARM::CHARM_HEADERS # charm++.h included in public header file
  PUBLIC error parallel performance
)


addCelloLib(parameters "parse.tab.c;lex.yy.c")
target_link_libraries(parameters
  PUBLIC error performance monitor parallel
  PRIVATE cello_component # comes from source files
)


addCelloLib(performance "")
target_link_libraries(performance
  PUBLIC error parallel parameters memory
)
if (use_papi)
  # it shouldn't take much to make this a PRIVATE dependency
  target_link_libraries(performance PUBLIC PAPI::papi)
endif()


addCelloLib(problem "")
target_link_libraries(problem
  PUBLIC cello_component charm_component error performance monitor memory
         parallel disk io parameters compute control mesh view data simulation
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)

addCelloLib(simulation "")
# somehow the main component needs to be a private dependency
target_link_libraries(simulation
  PUBLIC cello_component charm_component error memory performance monitor
         parallel disk problem mesh view data parameters io compute
  PRIVATE cello_component_PCH # (b/c cello.hpp included in ALL source files)
)


addCelloLib(test_Unit "")
target_link_libraries(test_Unit PRIVATE cello_component error parallel monitor)

#===============================================
# Defining Cello-Component Targets step 3 (of 3)
#===============================================
# define "main-targets" that define the necessary function to drive unit tests
# (each unit test needs to be linked against one of these targets). Hereafter,
# we refer to these targets as "tester"
#
# The different tester targets all define the main chare, but include different
# subsets of chare declarations
#
# TODO: we can DEFINITELY avoid adding some charm module dependencies to
# these targets

function(makeTester TARGET SRC)
  # construct a target that is linked into unit tests
  add_library(${TARGET} STATIC ${SRC})
  target_link_libraries(${TARGET}
    PUBLIC parallel cello_component monitor
    PRIVATE Boost::filesystem
    INTERFACE test_Unit # the source file compiled in this target doesn't
                        # directly depend on this, but all tests depend on this
  )
  setupCelloLibHelper(${TARGET})

endfunction(makeTester)

# this does not contain any chare declarations
makeTester(tester_default "main.cpp")


# this includes the chare declarations for the Mesh chares
makeTester(tester_mesh "main_mesh.cpp")
target_link_libraries(tester_mesh
  PRIVATE simulation # once we stop providing dummy implementations of methods
                     # for Block and Simulation in main.cpp, we can drop this
                     # dependency
)


# this includes the declarations for the mesh and simulation chares
makeTester(tester_simulation "main_simulation.cpp")
target_link_libraries(tester_simulation
  PRIVATE simulation # once we stop providing dummy implementations of methods
                     # for Block and Simulation in main.cpp, we can drop this
                     # dependency
)

# NOTE: main_enzo.cpp includes declarations for the mesh and simulation chares
# as well as chares from the Enzo-layer. This is only touched in the Enzo layer
# and is used to drive the actual simulation

# TODO: Consider making a target called something like CELLO_COMBINED that
# could be used in the enzo-layer (rather than linking against every
# individual component)

#==========================
# Define Unit-Test Binaries
#==========================

function(addUnitTestBinary BINNAME SRCS LIB TESTER)
  # Arguments:
  #   BINNAME:   Name of the test binary to build
  #   SRCS:      Source files
  #   LIB:       Library required
  #   TESTER:    Tester target library (usually, this is tester_default,
  #              tester_mesh, OR tester_simulation)

  if (BUILD_TESTING)
    # only declare the binary for a unit-test when BUILD_TESTING is ON (the
    # default case)
    add_executable(${BINNAME} ${SRCS})
    target_link_libraries(${BINNAME} PUBLIC ${LIB} ${TESTER})
    target_link_options(${BINNAME} PRIVATE ${Cello_TARGET_LINK_OPTIONS})
  endif()
endfunction(addUnitTestBinary)

# tests of the cello-component
addUnitTestBinary(test_type "test_Type.cpp" cello_component tester_default)

# tests of the disk component
addUnitTestBinary(test_disk_utils "test_DiskUtils.cpp" disk tester_default)
addUnitTestBinary(test_file_hdf5 "test_FileHdf5.cpp" disk tester_default)

# tests of the error component
addUnitTestBinary(test_error "test_Error.cpp" error tester_default)

# tests of the data component
addUnitTestBinary(test_particle "test_Particle.cpp" data tester_default)
addUnitTestBinary(test_scalar "test_Scalar.cpp" data tester_default)
addUnitTestBinary(test_field_data "test_FieldData.cpp" data tester_default)
addUnitTestBinary(test_field_descr "test_FieldDescr.cpp" data tester_default)
addUnitTestBinary(test_field "test_Field.cpp" data tester_default)
addUnitTestBinary(test_field_face "test_FieldFace.cpp" data tester_simulation)
addUnitTestBinary(test_grouping "test_Grouping.cpp" data tester_default)
addUnitTestBinary(test_itindex "test_ItIndex.cpp" data tester_simulation)

# tests of the io component
addUnitTestBinary(test_colormap "test_Colormap.cpp" io tester_simulation)
addUnitTestBinary(test_schedule "test_Schedule.cpp" io tester_simulation)

# assorted tests - we should probably divide this up by component
addUnitTestBinary(test_data "test_Data.cpp" mesh tester_mesh)
addUnitTestBinary(test_it_face "test_ItFace.cpp" mesh tester_mesh)
addUnitTestBinary(test_it_child "test_ItChild.cpp" mesh tester_mesh)
addUnitTestBinary(test_block_trace "test_BlockTrace.cpp" mesh tester_mesh)
addUnitTestBinary(test_sync "test_Sync.cpp" mesh tester_mesh)
addUnitTestBinary(test_face "test_Face.cpp" mesh tester_mesh)
addUnitTestBinary(test_face_fluxes "test_FaceFluxes.cpp" mesh tester_mesh)
addUnitTestBinary(test_flux_data "test_FluxData.cpp" mesh tester_mesh)
addUnitTestBinary(test_prolong_linear "test_ProlongLinear.cpp" mesh tester_mesh)
addUnitTestBinary(test_refresh "test_Refresh.cpp" mesh tester_mesh)
addUnitTestBinary(test_mask "test_Mask.cpp" mesh tester_mesh)
addUnitTestBinary(test_value "test_Value.cpp" mesh tester_mesh)
addUnitTestBinary(test_box "test_Box.cpp" mesh tester_mesh)
addUnitTestBinary(test_adapt "test_Adapt.cpp" mesh tester_mesh)

# test of the memory component
addUnitTestBinary(test_memory "test_Memory.cpp" memory tester_default)

# test of the monitor component
addUnitTestBinary(test_monitor "test_Monitor.cpp" monitor tester_default)

# tests of the parameters component
addUnitTestBinary(
  test_parameters "test_Parameters.cpp" parameters tester_default
)
addUnitTestBinary(test_parse "test_Parse.cpp" parameters tester_default)

# tests of the performance component
addUnitTestBinary(
  test_performance "test_Performance.cpp" performance tester_default
)
addUnitTestBinary(test_timer "test_Timer.cpp" performance tester_default)
if (use_papi)
  addUnitTestBinary(test_papi "test_Papi.cpp" performance tester_default)
endif()

# tests of the view component
addUnitTestBinary(test_view_collec "test_ViewCollec.cpp" view tester_default)
addUnitTestBinary(test_cello_view "test_CelloView.cpp" view tester_default)
addUnitTestBinary(
  test_string_ind_rd_only_map "test_StringIndRdOnlyMap.cpp" view tester_default
)

# the following test is broken (an API change is not reflected in the test)
#addUnitTestBinary(test_index "test_Index.cpp" mesh tester_mesh)

