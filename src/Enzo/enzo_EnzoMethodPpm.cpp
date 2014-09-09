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

void EnzoMethodPpm::compute ( CommBlock * comm_block) throw()
{
  if (!comm_block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  enzo_block->SolveHydroEquations 
    ( comm_block->cycle(), comm_block->time(), comm_block->dt() );
}

//----------------------------------------------------------------------

double EnzoMethodPpm::timestep ( CommBlock * comm_block ) throw()
{

  EnzoBlock * enzo_comm_block = static_cast<EnzoBlock*> (comm_block);
  Field field = enzo_comm_block->block()->field();

  enzo_float * density    = (enzo_float *)field.values("density");
  enzo_float * velocity_x = (enzo_float *)field.values("velocity_x");
  enzo_float * velocity_y = (enzo_float *)field.values("velocity_y");
  enzo_float * velocity_z = (enzo_float *)field.values("velocity_z");

  enzo_float a = 1, dadt;
  
  if (EnzoBlock::ComovingCoordinates)
    enzo_comm_block->CosmologyComputeExpansionFactor
      (enzo_comm_block->Time(), &a, &dadt);

  enzo_float dtBaryons      = ENZO_HUGE_VAL;
  enzo_float dtExpansion    = ENZO_HUGE_VAL;

  /* Compute the pressure. */

  EnzoMethodPressure method_pressure (EnzoBlock::Gamma);
  method_pressure.compute(enzo_comm_block);

  enzo_float * pressure = (enzo_float *) field.values("pressure");

 
  /* 2) Calculate dt from particles. */
 
  /* 3) Find dt from expansion. */
 
  if (EnzoBlock::ComovingCoordinates)
    if (enzo_comm_block->CosmologyComputeExpansionTimestep(enzo_comm_block->Time(), &dtExpansion) == ENZO_FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(ENZO_FAIL);
    }
 
  //   /* 4) Calculate minimum dt due to acceleration field (if present). */
 
  /* 5) calculate minimum timestep */

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
