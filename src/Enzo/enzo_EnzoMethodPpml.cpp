// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpml class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpml::EnzoMethodPpml( Parameters * parameters )
  : Hyperbolic(parameters)
{
}

//----------------------------------------------------------------------

void EnzoMethodPpml::compute_block
(
 Block * block,
 double t,
 double dt
 ) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveMHDEquations ( enzo_block->CycleNumber, dt);
}

//----------------------------------------------------------------------
