// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/cosmology/cosmology.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for the cosmology subcomponent of the \ref Enzo layer

#ifndef ENZO_COSMOLOGY_COSMOLOGY_HPP
#define ENZO_COSMOLOGY_COSMOLOGY_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Method, Physics

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock, EnzoEFltArrayMap

//----------------------------------------------------------------------
// Assorted Public Functions
//----------------------------------------------------------------------

// TODO: potentially move this into its own private header file (or just
//       declare it in the source file where its used)
extern "C" void FORTRAN_NAME(expand_terms)(
   int *rank, int *isize, int *idual, enzo_float *coef,
   int *imethod, enzo_float *gamma,
   enzo_float *p,  enzo_float *d, enzo_float *e, enzo_float *ge,
   enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dold, enzo_float *eold, enzo_float *geold,
   enzo_float *uold, enzo_float *vold, enzo_float *wold,
   int *icr, enzo_float *ecr, enzo_float *ecrold);

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "cosmology/EnzoPhysicsCosmology.hpp"
#include "cosmology/EnzoMethodComovingExpansion.hpp"
#include "cosmology/EnzoMethodCosmology.hpp"

#endif /* ENZO_COSMOLOGY_COSMOLOGY_HPP */
