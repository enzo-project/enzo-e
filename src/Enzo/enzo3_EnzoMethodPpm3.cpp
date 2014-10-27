// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm3.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm3 class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpm3::EnzoMethodPpm3 () : Method()
{
  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm3::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodPpm3::compute ( CommBlock * comm_block) throw()
{
  if (!comm_block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  enzo_block->SolveHydroEquations3
    ( comm_block->time(), comm_block->dt() );

}

//----------------------------------------------------------------------

double EnzoMethodPpm3::timestep ( CommBlock * comm_block ) throw()
{
  double dt = std::numeric_limits<double>::max();

  INCOMPLETE("EnzoMethodPpm3::timestep()");

//       FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
// 			     GridDimension+2,
// 			     GridStartIndex, GridEndIndex,
// 			     GridStartIndex+1, GridEndIndex+1,
// 			     GridStartIndex+2, GridEndIndex+2,
// 			     &HydroMethod, &ZEUSQuadraticArtificialViscosity,
// 			     CellWidth[0], CellWidth[1], CellWidth[2],
// 			     GridVelocity, GridVelocity+1, GridVelocity+2,
// 			     &Gamma, &PressureFree, &afloat,
// 			     BaryonField[DensNum], pressure_field,
// 			     BaryonField[Vel1Num], BaryonField[Vel2Num],
// 			     BaryonField[Vel3Num], &dtBaryons, &dtViscous);


  return dt;
}

