
include(charm)
addCharmModule( "enzo" )

file(GLOB SRC_FILES enzo*.cpp *F *enzo*hpp fortran.h fortran_types.h)

# TODO(pgrete) Check if that file can be removed for good
# -> PPML_sub.F:(.text+0x5ebf): multiple definition of `qdd6_'; CMakeFiles/enzo-e.dir/PPML_QDD6.F.o:PPML_QDD6.F:(.text+0x0): first defined here
#/usr/bin/ld: CMakeFiles/enzo-e.dir/enzo_EnzoMethodStarMakerStochasticSF.cpp.o: in function `EnzoMethodStarMakerStochasticSF::EnzoMethodStarMakerStochasticSF()':
list(FILTER SRC_FILES EXCLUDE REGEX "PPML_QDD6")

# TODO(pgrete) not sure about the current state of the src/Enzo/enzo_EnzoMethodStarMakerFIRE2.cpp file
# Currently is also contains the same names as StarMakerSF so there exists a conflict.
# Exluding it for now.
list(FILTER SRC_FILES EXCLUDE REGEX "FIRE2")

# TODO(pgrete) enzo_EnzoMethodGrackle.cpp:(.text+0x1f5): undefined reference to `vtable for EnzoMethodGrackle'
# so probably need to include higher level logic (for runs without Grackle)
list(FILTER SRC_FILES EXCLUDE REGEX "Grackle")


add_executable(enzo-e ${CMAKE_CURRENT_SOURCE_DIR}/../Cello/main_enzo.cpp ${SRC_FILES})
add_dependencies(enzo-e enzoCharmModule main_enzoCharmModule simulationCharmModule)
target_link_libraries(enzo-e PUBLIC charm simulation problem compute control mesh data io disk memory parameters error monitor parallel performance test_Unit cello pngwriter Boost::filesystem PNG::PNG HDF5::HDF5)
# TOOD(pgrete) need a cleaner place for the generated headers (current in Cello folder) or export header location
target_include_directories (enzo-e PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR}/../Cello ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
