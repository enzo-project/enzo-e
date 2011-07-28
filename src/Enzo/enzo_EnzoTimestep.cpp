// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestep.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoTimestep class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoTimestep::EnzoTimestep () throw()
  : Timestep()
{

}

//----------------------------------------------------------------------

double EnzoTimestep::compute ( const FieldDescr * field_descr,
			       Block * block ) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  FieldBlock * field_block = enzo_block->field_block();
  enzo_float * density_field    = 0;
  enzo_float * velocity_x_field = 0;
  enzo_float * velocity_y_field = 0;
  enzo_float * velocity_z_field = 0;

  density_field = (enzo_float *)field_block->field_values(enzo::field_density);
  if (enzo::GridRank >= 1)   velocity_x_field = 
    (enzo_float *)field_block->field_values(enzo::field_velocity_x);
  if (enzo::GridRank >= 2) velocity_y_field = 
    (enzo_float *)field_block->field_values(enzo::field_velocity_y);
  if (enzo::GridRank >= 3) velocity_z_field = 
    (enzo_float *)field_block->field_values(enzo::field_velocity_z);

  enzo_float a = 1, dadt;
  if (enzo::ComovingCoordinates)
    enzo_block->CosmologyComputeExpansionFactor(enzo_block->Time, &a, &dadt);
  //  enzo_float dt, dtTemp;
  enzo_float dtBaryons      = HUGE_VALF;
  //  enzo_float dtViscous      = HUGE_VALF;
  //  enzo_float dtParticles    = HUGE_VALF;
  enzo_float dtExpansion    = HUGE_VALF;
  //  enzo_float dtAcceleration = HUGE_VALF;


  /* Compute the pressure. */

  int nx,ny,nz;
  field_block -> size (&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr->ghosts(enzo::field_density,&gx,&gy,&gz);
  int mx,my,mz;
  mx = nx + 2*gx;
  my = ny + 2*gx;
  mz = nz + 2*gz;

  int size = mx*my*mz;

  // @@@ WARNING: should have temporary Field capability

  enzo_float * pressure_field = new enzo_float[size];
  for (int i=0; i<size; i++) pressure_field[i] = 0;

  int  result;
  if (enzo::DualEnergyFormalism)
    result = enzo_block->ComputePressureDualEnergyFormalism
      (enzo_block->Time, pressure_field);
  else
    result = enzo_block->ComputePressure
      (enzo_block->Time, pressure_field);
 
  if (result == ENZO_FAIL) {
    fprintf(stderr, "Error in grid->ComputePressure.\n");
    exit(EXIT_FAILURE);
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
 
  if (enzo::ComovingCoordinates)
    if (enzo_block->CosmologyComputeExpansionTimestep(enzo_block->Time, &dtExpansion) == ENZO_FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(ENZO_FAIL);
    }
 
  //   /* 4) Calculate minimum dt due to acceleration field (if present). */
 
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

  FORTRAN_NAME(calc_dt)(&enzo::GridRank, enzo_block->GridDimension, enzo_block->GridDimension+1,
			enzo_block->GridDimension+2,
			enzo_block->GridStartIndex, enzo_block->GridEndIndex,
			enzo_block->GridStartIndex+1, enzo_block->GridEndIndex+1,
			enzo_block->GridStartIndex+2, enzo_block->GridEndIndex+2,
			&enzo_block->CellWidth[0], 
			&enzo_block->CellWidth[1], 
			&enzo_block->CellWidth[2],
			&enzo::Gamma, &enzo::PressureFree, &a,
			density_field, pressure_field,
			velocity_x_field, 
			velocity_y_field, 
			velocity_z_field, 
			&dtBaryons);

  dtBaryons *= enzo::CourantSafetyNumber;
 
  double dt = dtBaryons;
  //  dt = MIN(dtBaryons, dtParticles);
  //  dt = MIN(dt, dtViscous);
  //   dt = MIN(dt, dtAcceleration);
  //  dt = MIN(dt, dtExpansion);

  delete [] pressure_field;

  return dt;
}

