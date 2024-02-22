// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/io/io.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for fluid-props subcomponent of the \ref Enzo layer

#ifndef ENZO_CHEMISTRY_CHEMISTRY_HPP
#define ENZO_CHEMISTRY_CHEMISTRY_HPP

// TODO: revisit this file in the future. We want to isolate <grackle.h> as a
//       private dependency (that isn't included in this public header)

#ifdef CONFIG_USE_GRACKLE
#include <stdlib.h>
extern "C" {
  #define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
  #include <grackle.h>
}
#else
extern "C" { // declare the names of Grackle types so can reduce the usage of
             // ifdef statements
  struct chemistry_data;
  struct chemistry_data_storage;
  struct code_units;
  struct grackle_field_data;
}
#endif

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/compute.hpp" // Compute
#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Method
#include "Cello/view.hpp" // CelloView

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock, EnzoEFltArrayMap

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "chemistry/GrackleChemistryData.hpp"
#include "chemistry/GrackleFacade.hpp"
#include "chemistry/EnzoComputeCoolingTime.hpp"
#include "chemistry/EnzoMethodGrackle.hpp"

#endif /* ENZO_CHEMISTRY_CHEMISTRY_HPP */
