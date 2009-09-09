/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    dt   - timestep
/
************************************************************************/
 
// Compute the timestep from all the constrains for this grid.
//
// Somebody fix the error handling in this routine! please.
//

#include "cello_hydro.h" 
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include "macros_and_parameters.h"
// #include "typedefs.h"
// #include "global_data.h"
// #include "Fluxes.h"
// #include "GridList.h"
// #include "ExternalBoundary.h"
// #include "Grid.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, 
		  float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);
 
 
float ComputeTimeStep()
{
 
  /* initialize */
 
  float dt, dtTemp;
  float dtBaryons      = HUGE_VAL;
  float dtViscous      = HUGE_VAL;
  float dtParticles    = HUGE_VAL;
  float dtExpansion    = HUGE_VAL;
  float dtAcceleration = HUGE_VAL;
  int dim, i, result;
 
  /* Compute the field size. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  float afloat = float(a);

  /* 1) Compute Courant condition for baryons. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }
 
    /* Compute the pressure. */
 
    float *pressure_field = new float[size];
    if (DualEnergyFormalism)
      result = ComputePressureDualEnergyFormalism(Time, pressure_field);
    else
      result = ComputePressure(Time, pressure_field);
 
    if (result == FAIL) {
      fprintf(stderr, "Error in grid->ComputePressure.\n");
      exit(EXIT_FAILURE);
    }
 
#ifdef UNUSED
    int Zero[3] = {0,0,0}, TempInt[3] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */
 
    /* Call fortran routine to do calculation. */
 
    FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
                               GridDimension+2,
//                        Zero, TempInt, Zero+1, TempInt+1, Zero+2, TempInt+2,
                          GridStartIndex, GridEndIndex,
                               GridStartIndex+1, GridEndIndex+1,
                               GridStartIndex+2, GridEndIndex+2,
			       CellWidth[0], CellWidth[1], CellWidth[2],
                               &Gamma, &PressureFree, &afloat,
                          BaryonField[DensNum], pressure_field,
                               BaryonField[Vel1Num], BaryonField[Vel2Num],
                               BaryonField[Vel3Num], &dtBaryons, &dtViscous);
 
    /* Clean up */
 
    delete pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= CourantSafetyNumber;
 
  }
 
  /* 2) Calculate dt from particles. */
 
//   if (NumberOfParticles > 0) {
 
//     /* Compute dt constraint from particle velocities. */
 
//     for (dim = 0; dim < GridRank; dim++) {
//       float dCell = CellWidth[dim][0]*a;
//       for (i = 0; i < NumberOfParticles; i++) {
//         dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
// 	dtParticles = min(dtParticles, dtTemp);
//       }
//     }
 
//     /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
//     dtParticles *= ParticleCourantSafetyNumber;
 
//   }
 
  /* 3) Find dt from expansion. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionTimestep(Time, &dtExpansion) == FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(FAIL);
    }
 
//   /* 4) Calculate minimum dt due to acceleration field (if present). */
 
//   if (SelfGravity) {
//     for (dim = 0; dim < GridRank; dim++)
//       if (AccelerationField[dim])
// 	for (i = 0; i < size; i++) {
// 	  dtTemp = sqrt(CellWidth[dim][0]/
// 			fabs(AccelerationField[dim][i])+tiny_number);
// 	  dtAcceleration = min(dtAcceleration, dtTemp);
// 	}
//     if (dtAcceleration != HUGE_VAL)
//       dtAcceleration *= 0.5;
//   }
 
  /* 5) calculate minimum timestep */
 
  dt = dtBaryons;
  //  dt = min(dtBaryons, dtParticles);
  //  dt = min(dt, dtViscous);
//   dt = min(dt, dtAcceleration);
//  dt = min(dt, dtExpansion);

 
  return dt;
}
