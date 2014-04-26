// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm () : Method()
{
  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodPpm::compute
(
 FieldDescr * field_descr, CommBlock * comm_block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);
  enzo_block->SolveHydroEquations 
    ( comm_block->cycle(), comm_block->time(), comm_block->dt() );
}

//----------------------------------------------------------------------

double EnzoMethodPpm::timestep
(
 const FieldDescr * field_descr,
 CommBlock *        comm_block
 ) const throw()
{

  EnzoBlock * enzo_comm_block = static_cast<EnzoBlock*> (comm_block);
  FieldBlock * field_block = enzo_comm_block->block()->field_block();
  enzo_float * density_field    = 0;
  enzo_float * velocity_x_field = 0;
  enzo_float * velocity_y_field = 0;
  enzo_float * velocity_z_field = 0;

  // DEBUG
    // double lower[3] = {0,0,0};
    // double upper[3] = {1,1,1};
    // enzo_comm_block->field_block()->print (field_descr,"dump",lower,upper);

  int index = enzo_comm_block->index(field_density);
  density_field = (enzo_float *)field_block->field_values(index);

  if (EnzoBlock::GridRank >= 1) {
    index = enzo_comm_block->index(field_velocity_x);
    velocity_x_field = (enzo_float *) field_block->field_values(index);
  }

  if (EnzoBlock::GridRank >= 2) {
    index = enzo_comm_block->index(field_velocity_y);
    velocity_y_field = (enzo_float *) field_block->field_values(index);
  }
  if (EnzoBlock::GridRank >= 3){
    index = enzo_comm_block->index(field_velocity_z);
    velocity_z_field = (enzo_float *) field_block->field_values(index);
  }

  enzo_float a = 1, dadt;
  
  if (EnzoBlock::ComovingCoordinates)
    enzo_comm_block->CosmologyComputeExpansionFactor(enzo_comm_block->Time(), &a, &dadt);
  //  enzo_float dt, dtTemp;
  enzo_float dtBaryons      = ENZO_HUGE_VAL;
  //  enzo_float dtViscous      = ENZO_HUGE_VAL;
  //  enzo_float dtParticles    = ENZO_HUGE_VAL;
  enzo_float dtExpansion    = ENZO_HUGE_VAL;
  //  enzo_float dtAcceleration = ENZO_HUGE_VAL;


  /* Compute the pressure. */

  int nx,ny,nz;
  field_block -> size (&nx,&ny,&nz);
  int gx,gy,gz;
  index = enzo_comm_block->index(field_density);
  field_descr->ghosts(index,&gx,&gy,&gz);
  int mx,my,mz;
  mx = nx + 2*gx;
  my = ny + 2*gx;
  mz = nz + 2*gz;

  int size = mx*my*mz;

  // @@@ WARNING: should have temporary Field capability

  enzo_float * pressure_field = new enzo_float[size];
  for (int i=0; i<size; i++) pressure_field[i] = 0;

  int  result;
  if (EnzoBlock::DualEnergyFormalism) {
    result = enzo_comm_block->ComputePressureDualEnergyFormalism
      (enzo_comm_block->Time(), pressure_field);
    //    printf ("%s:%d pressure %f\n",__FILE__,__LINE__,pressure_field[0]);
  }
  else {
    result = enzo_comm_block->ComputePressure
      (enzo_comm_block->Time(), pressure_field);
    //    printf ("%s:%d pressure %f\n",__FILE__,__LINE__,pressure_field[0]);
  }
 
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
 
  if (EnzoBlock::ComovingCoordinates)
    if (enzo_comm_block->CosmologyComputeExpansionTimestep(enzo_comm_block->Time(), &dtExpansion) == ENZO_FAIL) {
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
  //     if (dtAcceleration != ENZO_HUGE_VAL)
  //       dtAcceleration *= 0.5;
  //   }
 
  /* 5) calculate minimum timestep */


  TRACE3("enzo_comm_block->GridStartIndex = %d %d %d",
	 enzo_comm_block->GridStartIndex[0],
	 enzo_comm_block->GridStartIndex[1],
	 enzo_comm_block->GridStartIndex[2]);

  FORTRAN_NAME(calc_dt)(&EnzoBlock::GridRank, enzo_comm_block->GridDimension, enzo_comm_block->GridDimension+1,
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
			&EnzoBlock::Gamma, &EnzoBlock::PressureFree, &a,
			density_field, pressure_field,
			velocity_x_field, 
			velocity_y_field, 
			velocity_z_field, 
			&dtBaryons);

  //  printf ("density_field = %f\n",density_field[0]);
  TRACE1 ("dtBaryons: %f",dtBaryons);
  //  printf("Courant = %f\n",EnzoBlock::CourantSafetyNumber);

  //  printf ("%s:%d %f\n",__FILE__,__LINE__,dtBaryons);

  dtBaryons *= EnzoBlock::CourantSafetyNumber;

  //  printf ("%s:%d %f\n",__FILE__,__LINE__,dtBaryons);
  
  double dt = dtBaryons;
  //  dt = MIN(dtBaryons, dtParticles);
  //  dt = MIN(dt, dtViscous);
  //   dt = MIN(dt, dtAcceleration);
  //  dt = MIN(dt, dtExpansion);

  delete [] pressure_field;

  if (dt<=0.0) {
    FieldBlock * field_block = comm_block->block()->field_block();
    field_block->print(comm_block->name().c_str(),true);
  }

  //  printf ("%s:%d %f\n",__FILE__,__LINE__,dt);
  return dt;
}
