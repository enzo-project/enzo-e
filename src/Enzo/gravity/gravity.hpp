// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/gravity/gravity.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for gravity subcomponent within the \ref Enzo layer

#ifndef ENZO_GRAVITY_GRAVITY_HPP
#define ENZO_GRAVITY_GRAVITY_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/compute.hpp" // Matrix, Solver
#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Method

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "gravity/EnzoPhysicsGravity.hpp"
#include "gravity/EnzoMethodBackgroundAcceleration.hpp"
#include "gravity/EnzoMethodGravity.hpp"
#include "gravity/EnzoMethodPmDeposit.hpp"

#include "gravity/matrix/EnzoMatrixDiagonal.hpp"
#include "gravity/matrix/EnzoMatrixIdentity.hpp"
#include "gravity/matrix/EnzoMatrixLaplace.hpp"

#include "gravity/EnzoComputeAcceleration.hpp"

#include "gravity/solvers/EnzoSolverBiCgStab.hpp"
#include "gravity/solvers/EnzoSolverCg.hpp"
#include "gravity/solvers/EnzoSolverDd.hpp"
#include "gravity/solvers/EnzoSolverDiagonal.hpp"
#include "gravity/solvers/EnzoSolverJacobi.hpp"
#include "gravity/solvers/EnzoSolverMg0.hpp"

#endif /* ENZO_GRAVITY_GRAVITY_HPP */
