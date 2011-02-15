// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      Grid_SetMinimumSupport.cpp
/// @author    Greg Bryan
/// @date      November, 1998
/// @brief     Set the energy to provide minimal pressure support

#include "cello.hpp"

#include "enzo.hpp"
 
/* function prototypes */
 
// int CosmologyComputeExpansionFactor(ENZO_FLOAT time, ENZO_FLOAT *a, ENZO_FLOAT *dadt);
 
 
int EnzoDescr::SetMinimumSupport(float &MinimumSupportEnergyCoefficient)
{
  if (NumberOfBaryonFields > 0) {
 
    const float pi = 3.14159;
 
    /* Compute cosmology factors. */
 
    ENZO_FLOAT a = 1, dadt;
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == ENZO_FAIL) {
	fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
	return ENZO_FAIL;
      }
    float CosmoFactor = 1.0/a;
 
    /* Determine the size of the grids. */
 
    int dim, size = 1, i;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Find the density, gas energy, velocities & total energy
       (where appropriate). */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == ENZO_FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return ENZO_FAIL;
    }
 
    /* Set minimum GE. */
 
    MinimumSupportEnergyCoefficient =
      GravitationalConstant/(4.0*pi) / (pi * (Gamma*(Gamma-1.0))) *
      CosmoFactor * MinimumPressureSupportParameter *
      CELLWIDTH(0,0) * CELLWIDTH(0,0);
 
 
    /* PPM: set GE. */
 
    if (DualEnergyFormalism == TRUE) {
      for (i = 0; i < size; i++)
	BaryonField[GENum][i] = MAX(BaryonField[GENum][i],
				    MinimumSupportEnergyCoefficient *
				    BaryonField[DensNum][i]);
      if (GridRank != 3) return ENZO_FAIL;
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] = 
	  MAX(BaryonField[GENum][i] + 0.5*
	      (BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i] +
	       BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i] +
	       BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i]),
	      BaryonField[TENum][i]);
								
    }
    else {
      fprintf(stderr, "not implemented.\n");
      return ENZO_FAIL;
    }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  return ENZO_SUCCESS;
}
