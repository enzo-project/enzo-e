// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestepPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoTimestepPpml class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoTimestepPpml::EnzoTimestepPpml () throw()
  : Timestep()
{

}

//----------------------------------------------------------------------

void EnzoTimestepPpml::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  Timestep::pup(p);
}

//----------------------------------------------------------------------

double EnzoTimestepPpml::evaluate ( const FieldDescr * field_descr,
				   CommBlock * comm_block ) throw()
{
 
  EnzoBlock * enzo_comm_block = static_cast<EnzoBlock*> (comm_block);

  /* initialize */
 
  enzo_float dt;
  // enzo_float dtTemp;
  enzo_float dtBaryons      = ENZO_HUGE_VAL;
  // enzo_float dtViscous      = ENZO_HUGE_VAL;
  // enzo_float dtParticles    = ENZO_HUGE_VAL;
  // enzo_float dtExpansion    = ENZO_HUGE_VAL;
  // enzo_float dtAcceleration = ENZO_HUGE_VAL;
  // int dim, i, result;
 
  /* Compute the field size. */
 
  // int size = 1;
  // for (dim = 0; dim < GridRank; dim++)
  //   size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  enzo_float a = 1, dadt;
  if (EnzoBlock::ComovingCoordinates)
    enzo_comm_block->CosmologyComputeExpansionFactor
      (enzo_comm_block->Time(), &a, &dadt);
  //  float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  if (EnzoBlock::NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    // float *pressure_field;

    // if (HydroMethod != PPML_Isothermal3D) {
    //   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
    // 					 Vel3Num, TENum) == FAIL) {
    // 	fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
    // 	exit(FAIL);
    //   }
 
    // /* Compute the pressure. */
 
    //   pressure_field = new float[size];
    //   if (DualEnergyFormalism)
    // 	result = this->ComputePressureDualEnergyFormalism(Time, pressure_field);
    //   else
    // 	result = this->ComputePressure(Time, pressure_field);
 
    //   if (result == FAIL) {
    // 	fprintf(stderr, "Error in grid->ComputePressure.\n");
    // 	exit(EXIT_FAILURE);
    //   }
    // }

    /* Call fortran routine to do calculation. */
 
    FieldBlock * field_block = enzo_comm_block->block()->field_block();
    enzo_float * d_field    = 0;
    enzo_float * vx_field = 0;
    enzo_float * vy_field = 0;
    enzo_float * vz_field = 0;
    enzo_float * bx_field = 0;
    enzo_float * by_field = 0;
    enzo_float * bz_field = 0;

    d_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_density));
    vx_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_velox));
    vy_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_veloy));
    vz_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_veloz));
    bx_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_bfieldx));
    by_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_bfieldy));
    bz_field = (enzo_float *)field_block->field_values(enzo_comm_block->index(field_bfieldz));


    // DEBUG
    // double lower[3] = {0,0,0};
    // double upper[3] = {1,1,1};
    // enzo_comm_block->field_block()->print (field_descr,"dump",lower,upper);

    FORTRAN_NAME(calc_dt_ppml)
      (enzo_comm_block->GridDimension, 
       enzo_comm_block->GridDimension+1, 
       enzo_comm_block->GridDimension+2,
       enzo_comm_block->GridStartIndex, 
       enzo_comm_block->GridEndIndex,
       enzo_comm_block->GridStartIndex+1, 
       enzo_comm_block->GridEndIndex+1,
       enzo_comm_block->GridStartIndex+2, 
       enzo_comm_block->GridEndIndex+2,
       &enzo_comm_block->CellWidth[0], 
       &enzo_comm_block->CellWidth[1], 
       &enzo_comm_block->CellWidth[2],
       d_field,
       vx_field, vy_field, vz_field,
       bx_field, by_field, bz_field,
       &dtBaryons);
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
 
    // if (HydroMethod != PPML_Isothermal3D)
    //   delete pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= EnzoBlock::CourantSafetyNumber;
    
  }
 
  /* 2) Calculate dt from particles. */
 
  // if (NumberOfParticles > 0) {
 
  //   /* Compute dt constraint from particle velocities. */
 
  //   for (dim = 0; dim < GridRank; dim++) {
  //     float dCell = CellWidth[dim][0]*a;
  //     for (i = 0; i < NumberOfParticles; i++) {
  //       dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
  // 	dtParticles = min(dtParticles, dtTemp);
  //     }
  //   }
 
  //   /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
  //   dtParticles *= ParticleCourantSafetyNumber;
 
  // }
 
  /* 3) Find dt from expansion. */
 
  // if (EnzoBlock::ComovingCoordinates)
  //   if (enzo_comm_block->CosmologyComputeExpansionTimestep
  // 	(enzo_comm_block->Time(), &dtExpansion) == ENZO_FAIL) {
  //     fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
  //     exit(EXIT_FAILURE);
  //   }
 
  /* 4) Calculate minimum dt due to acceleration field (if present). */
 
  // if (SelfGravity) {
  //   for (dim = 0; dim < GridRank; dim++)
  //     if (AccelerationField[dim])
  // 	for (i = 0; i < size; i++) {
  // 	  dtTemp = sqrt(CellWidth[dim][0]/
  // 			fabs(AccelerationField[dim][i])+tiny_number);
  // 	  dtAcceleration = min(dtAcceleration, dtTemp);
  // 	}
  //   if (dtAcceleration != huge_number)
  //     dtAcceleration *= 0.5;
  // }
 
  /* 5) calculate minimum timestep */

  dt = std::numeric_limits<enzo_float>::max();

  dt = MIN(dt, dtBaryons);

  // dt = MIN(dt, dtParticles);
  // dt = min(dt, dtViscous);
  // dt = min(dt, dtAcceleration);
  // dt = min(dt, dtExpansion);
 
  /* Debugging info. */
 
  // if (debug1) {
  //   printf("ComputeTimeStep = %"FSYM" (", dt);
  //   if (NumberOfBaryonFields > 0)
  //     printf("Bar = %"FSYM" ", dtBaryons);
  //   if (HydroMethod == Zeus_Hydro)
  //     printf("Vis = %"FSYM" ", dtViscous);
  //   if (ComovingCoordinates)
  //     printf("Exp = %"FSYM" ", dtExpansion);
  //   if (dtAcceleration != huge_number)
  //     printf("Acc = %"FSYM" ", dtAcceleration);
  //   if (NumberOfParticles)
  //     printf("Part = %"FSYM" ", dtParticles);
  //   printf(")\n");
  // }
 

  return dt;
}

