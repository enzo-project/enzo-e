// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/mesh/mesh.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for mesh subcomponent within the \ref Enzo component

#ifndef ENZO_MESH_MESH_HPP
#define ENZO_MESH_MESH_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"
#include "Cello/charm.hpp"

#include "Cello/mesh.hpp"    // Block, Refine
#include "Cello/problem.hpp" // Prolong, Restrict

#include "Enzo/enzo.hpp" // enzo_float

//----------------------------------------------------------------------
// Assorted Public Functions
//----------------------------------------------------------------------

// TODO: this could probably be put into a separate header (and probably
// doesn't need to be publicly exposed)
extern "C" void FORTRAN_NAME(interpolate)
  (int *rank, enzo_float *pfield, int pdim[],
   int pis[], int pie[], int r[],
   enzo_float *field, int dim[], int is[], enzo_float *work,
   int *imethod, int *posflag,
   int *ierror);

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Enzo/mesh/EnzoRefineShock.hpp"
#include "Enzo/mesh/EnzoRefineParticleMass.hpp"
#include "Enzo/mesh/EnzoRefineMass.hpp"

#include "Enzo/mesh/EnzoProlong.hpp"
#include "Enzo/mesh/EnzoProlongMC1.hpp"
#include "Enzo/mesh/EnzoProlongPoisson.hpp"
#include "Enzo/mesh/EnzoRestrict.hpp"


#endif /* ENZO_MESH_MESH_HPP */
