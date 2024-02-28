# See LICENSE_CELLO file for license and copyright information

# create Enzo::tests, which represents defines source files for initializers
# that explicitly exist for testing purposes
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_tests
  tests.hpp
  EnzoInitialAccretionTest.cpp EnzoInitialAccretionTest.hpp
  EnzoInitialBBTest.cpp EnzoInitialBBTest.hpp
  EnzoInitialFeedbackTest.cpp EnzoInitialFeedbackTest.hpp
  EnzoInitialGrackleTest.cpp EnzoInitialGrackleTest.hpp
  EnzoInitialMergeSinksTest.cpp EnzoInitialMergeSinksTest.hpp
  EnzoInitialPpmlTest.cpp EnzoInitialPpmlTest.hpp
)
add_library(Enzo::tests ALIAS Enzo_tests)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_tests PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_tests PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_tests PRIVATE ${Cello_TARGET_LINK_OPTIONS})

