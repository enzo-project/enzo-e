/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
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
 
#include "cello_hydro.h"
 
/* function prototypes */
 
void my_exit(int status);

int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
			     int *ihydro, float *C2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
                             float *vgz, float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);
 
extern "C" void FORTRAN_NAME(calc_dt_ppml)(
                  int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz,
                  float *dn, float *vx, float *vy, float *vz, 
                             float *bx, float *by, float *bz, 
			     float *dt);
 
 
float ComputeTimeStep()
{
 
  /* initialize */
 
  float dt;
//   float dtTemp;
  float dtBaryons      = HUGE_VAL;
  float dtViscous      = HUGE_VAL;
  float dtParticles    = HUGE_VAL;
  float dtExpansion    = HUGE_VAL;
  float dtAcceleration = HUGE_VAL;
  int dim;
//   int i, result;
 
  /* Compute the field size. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
//   float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
//     int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
//     float *pressure_field;

//     if (HydroMethod != PPML_Isothermal3D) {
//       if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
// 					 Vel3Num, TENum) == FAIL) {
// 	fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
// 	exit(FAIL);
//       }
 
//     /* Compute the pressure. */
 
//       pressure_field = new float[size];
//       if (DualEnergyFormalism)
// 	result = this->ComputePressureDualEnergyFormalism(Time, pressure_field);
//       else
// 	result = this->ComputePressure(Time, pressure_field);
 
//       if (result == FAIL) {
// 	fprintf(stderr, "Error in grid->ComputePressure.\n");
// 	exit(EXIT_FAILURE);
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
// 	my_exit(EXIT_FAILURE);
      FORTRAN_NAME(calc_dt_ppml)
	(GridDimension, GridDimension+1, GridDimension+2,
	 GridStartIndex, GridEndIndex,
	 GridStartIndex+1, GridEndIndex+1,
	 GridStartIndex+2, GridEndIndex+2,
	 CellWidth[0], CellWidth[1], CellWidth[2],
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
 
  /* 4) Calculate minimum dt due to acceleration field (if present). */
 
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
 
  dt = min(dtBaryons, dtParticles);
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);
 
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
