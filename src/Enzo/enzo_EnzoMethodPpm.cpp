// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm 
(
 const FieldDescr * field_descr,
 EnzoConfig * enzo_config
) 
  : Method(),
    comoving_coordinates_(enzo_config->physics_cosmology)
{
  // Initialize the Refresh object

  refresh_ = new Refresh (4,0);
  refresh_->add_all_fields(enzo_config->num_fields);
  refresh_->set_sync_type(sync_neighbor);
  // refresh_list_.resize(1);

  // refresh_list_[0].set_ghost_depth(4);
  // refresh_list_[0].set_min_face_rank(0);
  // refresh_list_[0].add_all_fields(field_descr->field_count());
  // refresh_list_[0].set_sync_type(sync_neighbor);

  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | comoving_coordinates_;
}

//----------------------------------------------------------------------

void EnzoMethodPpm::compute ( Block * block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  if (block->is_leaf()) {

    enzo_block->SolveHydroEquations 
      ( block->time(), block->dt(), comoving_coordinates_ );

  }

  block->compute_done(); 
  
}

//----------------------------------------------------------------------

double EnzoMethodPpm::timestep ( Block * block ) const throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  enzo_float a = 1, dadt;
  
  if (comoving_coordinates_)
    enzo_block->CosmologyComputeExpansionFactor
      (enzo_block->time(), &a, &dadt);

  enzo_float dtBaryons      = ENZO_HUGE_VAL;
  enzo_float dtExpansion    = ENZO_HUGE_VAL;

  /* Compute the pressure. */

  EnzoComputePressure compute_pressure (EnzoBlock::Gamma,
					comoving_coordinates_);
  compute_pressure.compute(enzo_block);

  Field field = enzo_block->data()->field();

  int rank = block->rank();

  enzo_float * density    = (enzo_float *)field.values("density");
  enzo_float * velocity_x = (rank >= 1) ? 
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ? 
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ? 
    (enzo_float *)field.values("velocity_z") : NULL;
  enzo_float * pressure = (enzo_float *) field.values("pressure");
 
  /* 2) Calculate dt from particles. */
 
  /* 3) Find dt from expansion. */
 
  if (comoving_coordinates_)
    if (enzo_block->CosmologyComputeExpansionTimestep(block->time(), &dtExpansion) == ENZO_FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(ENZO_FAIL);
    }
 
  //   /* 4) Calculate minimum dt due to acceleration field (if present). */
 
  /* 5) calculate minimum timestep */

  FORTRAN_NAME(calc_dt)(&rank, 
			enzo_block->GridDimension, 
			enzo_block->GridDimension+1,
			enzo_block->GridDimension+2,
			enzo_block->GridStartIndex, 
			enzo_block->GridEndIndex,
			enzo_block->GridStartIndex+1, 
			enzo_block->GridEndIndex+1,
			enzo_block->GridStartIndex+2, 
			enzo_block->GridEndIndex+2,
			&enzo_block->CellWidth[0], 
			&enzo_block->CellWidth[1], 
			&enzo_block->CellWidth[2],
			&EnzoBlock::Gamma, &EnzoBlock::PressureFree, &a,
			density, pressure,
			velocity_x, 
			velocity_y, 
			velocity_z, 
			&dtBaryons);

  TRACE1 ("dtBaryons: %f",dtBaryons);

  dtBaryons *= EnzoBlock::CourantSafetyNumber;

  double dt = dtBaryons;

  return dt;
}
