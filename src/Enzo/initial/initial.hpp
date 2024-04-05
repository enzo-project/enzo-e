// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/initial/initial.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for mesh subcomponent within the \ref Enzo component

#ifndef ENZO_INITIAL_INITIAL_HPP
#define ENZO_INITIAL_INITIAL_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"
#include "Cello/charm.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Initial

#include "Enzo/enzo.hpp" // enzo_float, EnzoConfig

//----------------------------------------------------------------------
// Assorted Public Functions
//----------------------------------------------------------------------

// TODO: potentially move these into their own private header file
extern "C" void FORTRAN_NAME(turboinit)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v, enzo_float *w,
   int *in, int *jn, int *kn,
   int *ig, int *jg, int *kg);

extern "C" void FORTRAN_NAME(turboinit2d)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v,
   int *in, int *jn,
   int *ig, int *jg);

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Enzo/initial/EnzoInitialBCenter.hpp"
#include "Enzo/initial/EnzoInitialCloud.hpp"
#include "Enzo/initial/EnzoInitialCollapse.hpp"
#include "Enzo/initial/EnzoInitialCosmology.hpp"
#include "Enzo/initial/EnzoInitialImplosion2.hpp"
#include "Enzo/initial/EnzoInitialInclinedWave.hpp"
#include "Enzo/initial/EnzoInitialSedovArray2.hpp"
#include "Enzo/initial/EnzoInitialSedovArray3.hpp"
#include "Enzo/initial/EnzoInitialSedovRandom.hpp"
#include "Enzo/initial/EnzoInitialShockTube.hpp"
#include "Enzo/initial/EnzoInitialSoup.hpp"
#include "Enzo/initial/EnzoInitialTurbulence.hpp"
#include "Enzo/initial/EnzoInitialIsolatedGalaxy.hpp"
#include "Enzo/initial/EnzoInitialBurkertBodenheimer.hpp"
#include "Enzo/initial/EnzoInitialShuCollapse.hpp"

#include "Enzo/initial/obsolete/EnzoInitialPm.hpp"

#endif /* ENZO_INITIAL_INITIAL_HPP */
