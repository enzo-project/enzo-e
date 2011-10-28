// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1: Alexei Kritsuk, May 2009 - added PPML/MHD timestep restriction
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
 
#include "cello.hpp"

#include "enzo.hpp"

enzo_float EnzoBlock::ComputeTimeStep()
{

  TRACE("");
  /* initialize */
 
  enzo_float dt;
//   enzo_float dtTemp;
  enzo_float dtBaryons      = HUGE_VAL;
  enzo_float dtViscous      = HUGE_VAL;
  enzo_float dtParticles    = HUGE_VAL;
  enzo_float dtExpansion    = HUGE_VAL;
  enzo_float dtAcceleration = HUGE_VAL;
  int dim;
//   int i, result;
 
  /* Compute the field size. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  enzo_float a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time(), &a, &dadt);
//   enzo_float afloat = enzo_float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
//     int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
//     enzo_float *pressure_field;

//     if (HydroMethod != PPML_Isothermal3D) {
//       if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
// 					 Vel3Num, TENum) == ENZO_FAIL) {
// 	fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
// 	exit(ENZO_FAIL);
//       }
 
//     /* Compute the pressure. */
 
//       pressure_field = new enzo_float[size];
//       if (DualEnergyFormalism)
// 	result = this->ComputePressureDualEnergyFormalism(Time(), pressure_field);
//       else
// 	result = this->ComputePressure(Time(), pressure_field);
 
//       if (result == ENZO_FAIL) {
// 	fprintf(stderr, "Error in grid->ComputePressure.\n");
// 	exit(EXIT_ENZO_FAILURE);
//       }
//     }

#ifdef UNUSED
    int Zero[3] = {0,0,0}, TempInt[3] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */
 
    /* Call fortran routine to do calculation. */
 
//     if (HydroMethod == PPML_Isothermal3D) {
//       if (GridRank !=3) 
// 	my_exit(EXIT_ENZO_FAILURE);
      FORTRAN_NAME(calc_dt_ppml)
	(GridDimension, GridDimension+1, GridDimension+2,
	 GridStartIndex, GridEndIndex,
	 GridStartIndex+1, GridEndIndex+1,
	 GridStartIndex+2, GridEndIndex+2,
	 &CellWidth[0], &CellWidth[1], &CellWidth[2],
	 BaryonField[0], 
	 BaryonField[1], BaryonField[2], BaryonField[3], 
	 BaryonField[4], BaryonField[5], BaryonField[6], 
	 &dtBaryons);

      //     }
 
//     else
//       FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
//                                GridDimension+2,
// //                        Zero, TempInt, Zero+1, TempInt+1, Zero+2, TempInt+2,
//                           GridStartIndex, GridEndIndex,
//                                GridStartIndex+1, GridEndIndex+1,
//                                GridStartIndex+2, GridEndIndex+2,
// 			       &HydroMethod, &ZEUSQuadraticArtificialViscosity,
//                           CellWidth[0], CellWidth[1], CellWidth[2],
//                                GridVelocity, GridVelocity+1, GridVelocity+2,
//                                &Gamma, &PressureFree, &afloat,
//                           BaryonField[DensNum], pressure_field,
//                                BaryonField[Vel1Num], BaryonField[Vel2Num],
//                                BaryonField[Vel3Num], &dtBaryons, &dtViscous);
 
    /* Clean up */
 
//     if (HydroMethod != PPML_Isothermal3D)
//       delete pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= CourantSafetyNumber;
    
  }
 
  /* 2) Calculate dt from particles. */
 
//   if (NumberOfParticles > 0) {
 
//     /* Compute dt constraint from particle velocities. */
 
//     for (dim = 0; dim < GridRank; dim++) {
//       enzo_float dCell = CellWidth[dim][0]*a;
//       for (i = 0; i < NumberOfParticles; i++) {
//         dtTemp = dCell/MAX(fabs(ParticleVelocity[dim][i]), tiny_number);
// 	dtParticles = MIN(dtParticles, dtTemp);
//       }
//     }
 
//     /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
//     dtParticles *= ParticleCourantSafetyNumber;
 
//   }
 
  /* 3) Find dt from expansion. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionTimestep(Time(), &dtExpansion) == ENZO_FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(ENZO_FAIL);
    }
 
  /* 4) Calculate minimum dt due to acceleration field (if present). */
 
//   if (SelfGravity) {
//     for (dim = 0; dim < GridRank; dim++)
//       if (AccelerationField[dim])
// 	for (i = 0; i < size; i++) {
// 	  dtTemp = sqrt(CellWidth[dim][0]/
// 			fabs(AccelerationField[dim][i])+tiny_number);
// 	  dtAcceleration = MIN(dtAcceleration, dtTemp);
// 	}
//     if (dtAcceleration != HUGE_VAL)
//       dtAcceleration *= 0.5;
//   }
 
  /* 5) calculate minimum timestep */
 
  dt = MIN(dtBaryons, dtParticles);
  dt = MIN(dt, dtViscous);
  dt = MIN(dt, dtAcceleration);
  dt = MIN(dt, dtExpansion);
 
  /* Debugging info. */
 
//   if (debug1) {
//     printf("ComputeTimeStep = %"FSYM" (", dt);
//     if (NumberOfBaryonFields > 0)
//       printf("Bar = %"FSYM" ", dtBaryons);
// //     if (HydroMethod == Zeus_Hydro)
// //       printf("Vis = %"FSYM" ", dtViscous);
//     if (ComovingCoordinates)
//       printf("Exp = %"FSYM" ", dtExpansion);
//     if (dtAcceleration != HUGE_VAL)
//       printf("Acc = %"FSYM" ", dtAcceleration);
//     if (NumberOfParticles)
//       printf("Part = %"FSYM" ", dtParticles);
//     printf(")\n");
//   }
 
  return dt;
}
