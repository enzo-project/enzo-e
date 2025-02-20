// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/chemistry/chemistry.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for chemistry subcomponent of the \ref Enzo layer

#ifndef ENZO_CHEMISTRY_CHEMISTRY_HPP
#define ENZO_CHEMISTRY_CHEMISTRY_HPP

// it's important to include this first
// - eventually, it would be great to finish refactoring so that the grackle
//   datatypes don't need to be declared in the headers
// - the placement of this include statement will be less important after we
//   stop including the this file in "_enzo.hpp"
#include "chemistry/chemistry_grackleincl.hpp"

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/compute.hpp"  // Compute
#include "Cello/mesh.hpp"     // Block
#include "Cello/problem.hpp"  // Method
#include "Cello/view.hpp"     // CelloView

#include "Enzo/enzo.hpp"  // enzo_float, EnzoBlock, EnzoEFltArrayMap

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "chemistry/GrackleChemistryData.hpp"
#include "chemistry/GrackleFacade.hpp"
#include "chemistry/EnzoComputeCoolingTime.hpp"
#include "chemistry/EnzoMethodGrackle.hpp"

#endif /* ENZO_CHEMISTRY_CHEMISTRY_HPP */
