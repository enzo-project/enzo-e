// $Id$
// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  COMPUTES THE MAXIMUM ALLOWED EXPANSION TIMESTEP AT GIVEN TIME
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include "cello.hpp"

#include "enzo.hpp"
 
//---------------------------------------------------------------------- 
 
int EnzoDescr::CosmologyComputeExpansionTimestep
(enzo_float time, enzo_float *dtExpansion)
{
 
  /* Error check. */
 
  if (InitialTimeInCodeUnits == 0) {
    fprintf(stderr, "The cosmology parameters seem to be improperly set.\n");
    return ENZO_FAIL;
  }
 
  /* Compute the expansion factors. */
 
  enzo_float a, dadt;
  if (CosmologyComputeExpansionFactor(time, &a, &dadt) == ENZO_FAIL) {
    fprintf(stderr, "Error in ComputeExpnasionFactors.\n");
    return ENZO_FAIL;
  }
 
  /* Compute the maximum allwed timestep given the maximum allowed
     expansion factor. */
 
  *dtExpansion = MaxExpansionRate*a/dadt;
 
  return ENZO_SUCCESS;
}
