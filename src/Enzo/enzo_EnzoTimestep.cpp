// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestep.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoTimestep class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoTimestep::EnzoTimestep (EnzoDescr * enzo) throw()
  : Timestep(),
    enzo_(enzo)
{

}

//----------------------------------------------------------------------

double EnzoTimestep::compute ( DataBlock * data_block ) throw()
{

  FieldBlock * field_block = data_block->field_block();
  float * density_field    = 0;
  float * velocity_x_field = 0;
  float * velocity_y_field = 0;
  float * velocity_z_field = 0;

  density_field    = (float *)field_block->field_values(enzo_->field_density);
  if (enzo_->GridRank >= 1)   velocity_x_field = 
    (float *)field_block->field_values(enzo_->field_velocity_x);
  if (enzo_->GridRank >= 2) velocity_y_field = 
    (float *)field_block->field_values(enzo_->field_velocity_y);
  if (enzo_->GridRank >= 3) velocity_z_field = 
    (float *)field_block->field_values(enzo_->field_velocity_z);

  ENZO_FLOAT a = 1, dadt;
  if (enzo_->ComovingCoordinates)
    enzo_->CosmologyComputeExpansionFactor(enzo_->Time, &a, &dadt);
  float afloat = float(a);
  //  float dt, dtTemp;
  float dtBaryons      = HUGE_VALF;
  float dtViscous      = HUGE_VALF;
  //  float dtParticles    = HUGE_VALF;
  float dtExpansion    = HUGE_VALF;
  //  float dtAcceleration = HUGE_VALF;


  /* Compute the pressure. */

  int nx,ny,nz;
  field_block -> size (&nx,&ny,&nz);
  int gx,gy,gz;
  field_block->field_descr()->ghosts(enzo_->field_density,&gx,&gy,&gz);
  int mx,my,mz;
  mx = nx + 2*gx;
  my = ny + 2*gx;
  mz = nz + 2*gz;

  int size = mx*my*mz;

  // @@@ WARNING: should have temporary Field capability

  float * pressure_field = new float[size];
  for (int i=0; i<size; i++) pressure_field[i] = 0;

  int  result;
  if (enzo_->DualEnergyFormalism)
    result = enzo_->ComputePressureDualEnergyFormalism(enzo_->Time, pressure_field);
  else
    result = enzo_->ComputePressure(enzo_->Time, pressure_field);
 
  if (result == ENZO_FAIL) {
    fprintf(stderr, "Error in grid->ComputePressure.\n");
    exit(EXIT_FAILURE);
  }
 
  
  /* 2) Calculate dt from particles. */
 
//   if (NumberOfParticles > 0) {
 
//     /* Compute dt constraint from particle velocities. */
 
//     for (dim = 0; dim < GridRank; dim++) {
//       float dCell = CellWidth[dim][0]*a;
//       for (i = 0; i < NumberOfParticles; i++) {
//         dtTemp = dCell/MAX(fabs(ParticleVelocity[dim][i]), tiny_number);
// 	dtParticles = MIN(dtParticles, dtTemp);
//       }
//     }
 
//     /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
//     dtParticles *= ParticleCourantSafetyNumber;
 
//   }
 
  /* 3) Find dt from expansion. */
 
  if (enzo_->ComovingCoordinates)
    if (enzo_->CosmologyComputeExpansionTimestep(enzo_->Time, &dtExpansion) == ENZO_FAIL) {
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

  printf("density = %p\n",density_field);
  FORTRAN_NAME(calc_dt)(&enzo_->GridRank, enzo_->GridDimension, enzo_->GridDimension+1,
			enzo_->GridDimension+2,
			enzo_->GridStartIndex, enzo_->GridEndIndex,
			enzo_->GridStartIndex+1, enzo_->GridEndIndex+1,
			enzo_->GridStartIndex+2, enzo_->GridEndIndex+2,
			enzo_->CellWidth[0], enzo_->CellWidth[1], enzo_->CellWidth[2],
			&enzo_->Gamma, &enzo_->PressureFree, &afloat,
			density_field, pressure_field,
			velocity_x_field, 
			velocity_y_field, 
			velocity_z_field, 
			&dtBaryons, &dtViscous);

  dtBaryons *= enzo_->CourantSafetyNumber;
 
  double dt = dtBaryons;
  //  dt = MIN(dtBaryons, dtParticles);
  //  dt = MIN(dt, dtViscous);
  //   dt = MIN(dt, dtAcceleration);
  //  dt = MIN(dt, dtExpansion);

  delete [] pressure_field;
  pressure_field = 0;

  return dt;
}

