// $Id: user_MethodEnzoTimestep.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoTimestep.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the MethodEnzoTimestep class

#include "enzo.hpp"
#include "data.hpp"
#include "user.hpp"
#include "parameters.hpp"

#include "cello_hydro.h"

//----------------------------------------------------------------------

MethodEnzoTimestep::MethodEnzoTimestep (Enzo * enzo) throw()
  : UserTimestep(),
    pressure_field_(0),
    afloat_(0),
    dtBaryons_(0),
    dtViscous_(0),
    dtExpansion_(0),
    enzo_(enzo)
{

}

//----------------------------------------------------------------------

void MethodEnzoTimestep::initialize (DataDescr * data_descr) throw()
{
}

//----------------------------------------------------------------------

void MethodEnzoTimestep::finalize ( DataDescr * data_descr ) throw ()
{
}

//----------------------------------------------------------------------

void MethodEnzoTimestep::initialize_block ( DataBlock * data_block ) throw ()
{
}

//----------------------------------------------------------------------

void MethodEnzoTimestep::finalize_block ( DataBlock * data_block ) throw ()
{
 
  delete [] pressure_field_;
  pressure_field_ = 0;
 
}

//----------------------------------------------------------------------

double MethodEnzoTimestep::compute_block ( DataBlock * data_block ) throw()
{

  FieldBlock * field_block = data_block->field_block();
  float * density_field    = (float *)field_block->field_values(enzo_->field_density);
  float * velocity_x_field = (float *)field_block->field_values(enzo_->field_velocity_x);
  float * velocity_y_field = (float *)field_block->field_values(enzo_->field_velocity_y);
  float * velocity_z_field = (float *)field_block->field_values(enzo_->field_velocity_z);

  ENZO_FLOAT a = 1, dadt;
  if (enzo_->ComovingCoordinates)
    enzo_->CosmologyComputeExpansionFactor(enzo_->Time, &a, &dadt);
  afloat_ = float(a);
  //  float dt, dtTemp;
  dtBaryons_      = HUGE_VALF;
  dtViscous_      = HUGE_VALF;
  //  float dtParticles    = HUGE_VALF;
  dtExpansion_    = HUGE_VALF;
  //  float dtAcceleration = HUGE_VALF;


  /* Compute the pressure. */

  int nx,ny,nz;
  field_block -> dimensions (&nx,&ny,&nz);
  int gx,gy,gz;
  field_block->field_descr()->ghosts(enzo_->field_density,&gx,&gy,&gz);
  int mx,my,mz;
  mx = nx + 2*gx;
  my = ny + 2*gx;
  mz = nz + 2*gz;

  int size = mx*my*mz;

  // @@@ WARNING: should have temporary Field capability
  // @@@ WARNING: this is dangerous if multiple blocks are being
  // @@@ WARNING: processed simultaneously, e.g. new block initialize
  // @@@ WARNING: is called before previous finalize.  ASSERT should
  // @@@ WARNING: prevent this from happening, but is overly restrictive

  ASSERT("MethodEnzoTimestep::initialize_block",
	 "new pressure field allocated before previous was deallocated",
	 pressure_field_ == NULL);
  pressure_field_ = new float[size];

  int  result;
  if (enzo_->DualEnergyFormalism)
    result = enzo_->ComputePressureDualEnergyFormalism(enzo_->Time, pressure_field_);
  else
    result = enzo_->ComputePressure(enzo_->Time, pressure_field_);
 
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
    if (enzo_->CosmologyComputeExpansionTimestep(enzo_->Time, &dtExpansion_) == ENZO_FAIL) {
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

  FORTRAN_NAME(calc_dt)(&enzo_->GridRank, enzo_->GridDimension, enzo_->GridDimension+1,
			enzo_->GridDimension+2,
			enzo_->GridStartIndex, enzo_->GridEndIndex,
			enzo_->GridStartIndex+1, enzo_->GridEndIndex+1,
			enzo_->GridStartIndex+2, enzo_->GridEndIndex+2,
			enzo_->CellWidth[0], enzo_->CellWidth[1], enzo_->CellWidth[2],
			&enzo_->Gamma, &enzo_->PressureFree, &afloat_,
			density_field, pressure_field_,
			velocity_x_field, 
			velocity_y_field, 
			velocity_z_field, 
			&dtBaryons_, &dtViscous_);

  dtBaryons_ *= enzo_->CourantSafetyNumber;
 
  double dt = dtBaryons_;
  //  dt = MIN(dtBaryons, dtParticles);
  //  dt = MIN(dt, dtViscous);
  //   dt = MIN(dt, dtAcceleration);
  //  dt = MIN(dt, dtExpansion);

 
  return dt;
}

