// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @todo     Remove unnecessary enzo_block-related parameters
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm ( Parameters * parameters )
  : Hyperbolic(parameters)
{
  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::compute_block
(
 Block * block,
 double t,
 double dt
 ) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveHydroEquations ( enzo_block->CycleNumber, dt);
}

