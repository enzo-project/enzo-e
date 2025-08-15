# See LICENSE_CELLO file for license and copyright information

# add source files related to the enzo target
#
# NOTE: this subdirectory was originally just named core, but that caused
#       charmc to produce an error message when it was called to translate a
#       charm++ interface file.
#       - The error message is benign. It occurs during cleanup in charmc when
#         it calls `rm -f core` (presumably to remove any core dump files), and
#         rm complains that it can't remove a directory.
#       - We choose to avoid that name in order to avoid worrying users


target_sources(enzo PRIVATE
  EnzoBlock.cpp
  EnzoBlock.hpp
  EnzoBoundary.cpp
  EnzoBoundary.hpp
  EnzoConfig.cpp
  EnzoConfig.hpp
  EnzoFactory.cpp
  EnzoFactory.hpp
  EnzoMethodBalance.cpp
  EnzoMethodBalance.hpp
  EnzoMsgCheck.cpp
  EnzoMsgCheck.hpp
  EnzoProblem.cpp
  EnzoProblem.hpp
  EnzoSimulation.cpp
  EnzoSimulation.hpp
  EnzoState.cpp
  EnzoState.hpp
  EnzoStopping.cpp
  EnzoStopping.hpp
  EnzoUnits.cpp
  EnzoUnits.hpp
  enzo.cpp
)

if (BUILD_TESTING)
  # Add a unit test
  add_executable(test_enzo_units test_EnzoUnits.cpp)
  target_link_libraries(test_enzo_units PRIVATE enzo main_enzo)
  target_link_options(test_enzo_units PRIVATE ${Cello_TARGET_LINK_OPTIONS})
endif()
