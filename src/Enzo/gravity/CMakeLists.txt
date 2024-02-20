# See LICENSE_CELLO file for license and copyright information

# create Enzo::gravity, which represents source files related to implementing
# gravity.
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_gravity
  gravity.hpp
  EnzoComputeAcceleration.cpp EnzoComputeAcceleration.hpp
  EnzoMethodBackgroundAcceleration.cpp EnzoMethodBackgroundAcceleration.hpp
  EnzoMethodGravity.cpp EnzoMethodGravity.hpp
  EnzoMethodPmDeposit.cpp EnzoMethodPmDeposit.hpp
  grid_cic.F

  matrix/EnzoMatrixDiagonal.cpp matrix/EnzoMatrixDiagonal.hpp
  matrix/EnzoMatrixIdentity.cpp matrix/EnzoMatrixIdentity.hpp
  matrix/EnzoMatrixLaplace.cpp matrix/EnzoMatrixLaplace.hpp

  solvers/EnzoSolverBiCgStab.cpp solvers/EnzoSolverBiCgStab.hpp
  solvers/EnzoSolverCg.cpp solvers/EnzoSolverCg.hpp
  solvers/EnzoSolverDd.cpp solvers/EnzoSolverDd.hpp
  solvers/EnzoSolverDiagonal.cpp solvers/EnzoSolverDiagonal.hpp
  solvers/EnzoSolverJacobi.cpp solvers/EnzoSolverJacobi.hpp
  solvers/EnzoSolverMg0.cpp solvers/EnzoSolverMg0.hpp
)
add_library(Enzo::gravity ALIAS Enzo_gravity)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_gravity PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_gravity PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_gravity PRIVATE ${Cello_TARGET_LINK_OPTIONS})
