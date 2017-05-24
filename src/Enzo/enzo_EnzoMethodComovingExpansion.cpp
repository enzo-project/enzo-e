// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodComovingExpansion.cpp
/// @author   Britton Smith (bds006@ucsd.edu)
/// @date     Wed May 24 12:25:56 PDT 2017
/// @brief    Implements comoving expansion class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodComovingExpansion::EnzoMethodComovingExpansion
(
 const FieldDescr * field_descr,
 EnzoConfig * enzo_config
) 
  : Method(),
    comoving_coordinates_(enzo_config->physics_cosmology)
{
  // do nothing!
}

//----------------------------------------------------------------------

void EnzoMethodComovingExpansion::compute ( Block * block) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  Field field = enzo_block->data()->field();

  /* Only do this if
     1. this is a leaf block
     2. we are using comoving coordinates
     3. baryon fields are present.
  */
  if !(block->is_leaf() &&
       comoving_coordinates_ &&
       field.field_count() > 0)
    { block->compute_done(); }

  /* Compute adot/a at time = t-1/2dt (time-centered). */

  enzo_float a, dadt, time, old_time;
  time     = enzo_block->time();
  old_time = enzo_block->old_time();

  if (enzo_block->CosmologyComputeExpansionFactor(
          0.5*(enzo_block->time()+enzo_block->old_time()),
          &a, &dadt) == ENZO_FAIL) {
    fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
    exit(ENZO_FAIL);
  }

  enzo_float Coefficient = enzo_block->dt()*dadt/a;

  /* Determine the size of the grids. */

  int i, dim, size = 1;
  int rank = enzo_block->rank();
  for (dim = 0; dim < rank; dim++)
    size *= GridDimension[dim];

  /* If we can, compute the pressure at the mid-point.
     We can, because we will always have an old baryon field now. */
  const int in = cello::index_static();
  enzo_float PressureTime = 0.5 * (time + old_time);
  enzo_float * pressure = new enzo_float[size];
  int rval;

  if (EnzoBlock::DualEnergyFormalism[in]) {
    rval = enzo_block->ComputePressureDualEnergyFormalism(
        PressureTime, pressure, comoving_coordinates_);
  }
  else{
    rval = enzo_block->ComputePressure(
        PressureTime, pressure, comoving_coordinates_);
  }
  if (rval == ENZO_FAIL) {
    fprintf(stderr, "Error in ComputePressureDualEnergyFormalism or ComputePressure.\n");
    exit(ENZO_FAIL);
  }

  // hard-code hydromethod for PPM for now
  int HydroMethod = 0;

  // hard-code CR method to off
  int CRModel = 0;
  enzo_float * cr_field, old_cr_field;
  cr_field = old_cr_field = NULL;

  /* Get the necessary fields. */
  // For now, assume the existence of "old_<field>".

  enzo_float * density         = (enzo_float *) field.values("density");
  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  enzo_float * internal_energy = (enzo_float *) field.values("internal_energy");
  enzo_float * old_density         = (enzo_float *) field.values("old_density");
  enzo_float * old_total_energy    = (enzo_float *) field.values("old_total_energy");
  enzo_float * old_internal_energy = (enzo_float *) field.values("old_internal_energy");

  enzo_float * velocity_x = NULL;
  enzo_float * velocity_y = NULL;
  enzo_float * velocity_z = NULL;
  enzo_float * old_velocity_x = NULL;
  enzo_float * old_velocity_y = NULL;
  enzo_float * old_velocity_z = NULL;

  velocity_x = (enzo_float *) field.values("velocity_x");
  old_velocity_x = (enzo_float *) field.values("old_velocity_x");
  if (rank >= 2) {
    velocity_y = (enzo_float *) field.values("velocity_y");
    old_velocity_y = (enzo_float *) field.values("old_velocity_y");
  }
  if (rank >= 3) {
    velocity_z = (enzo_float *) field.values("velocity_z");
    old_velocity_z = (enzo_float *) field.values("old_velocity_z");
  }

  /* Call fortran routine to do the real work. */

  FORTRAN_NAME(expand_terms)(
      &rank, &size, &EnzoBlock::DualEnergyFormalism[in], &Coefficient,
      (int*) &HydroMethod, &EnzoBlock::Gamma[in],
      pressure,
      density, total_energy, internal_energy,
      velocity_x, velocity_y, velocity_z,
      old_density, old_total_energy, old_internal_energy,
      old_velocity_x, old_velocity_y, old_velocity_z,
      &CRModel, cr_field, old_cr_field);

  block->compute_done();

}

//----------------------------------------------------------------------

double EnzoMethodComovingExpansion::timestep( Block * block ) const throw()
{

  enzo_float dtExpansion = ENZO_HUGE_VAL;

  if !(comoving_coordinates_)
    return (double) dtExpansion;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  enzo_float a = 1, dadt;

  if (enzo_block->CosmologyComputeExpansionFactor(
                  enzo_block->time(), &a, &dadt) == ENZO_FAIL) {
    fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
    exit(ENZO_FAIL);
  }

  if (enzo_block->CosmologyComputeExpansionTimestep(
                  enzo_block->time(), &dtExpansion) == ENZO_FAIL) {
    fprintf(stderr, "Error in ComputeExpansionTimestep.\n");
    exit(ENZO_FAIL);
  }

  return dtExpansion;

}

//----------------------------------------------------------------------

void EnzoMethodComovingExpansion::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

//======================================================================
