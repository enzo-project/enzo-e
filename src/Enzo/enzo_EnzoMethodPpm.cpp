// $Id$
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
 Parameters * parameters,
 EnzoDescr * enzo
 )
  : Hyperbolic(parameters),
    enzo_(enzo)
{
  // PPM parameters initialized in EnzoDescr::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPpm::compute_block
(
 DataBlock * data_block,
 double t,
 double dt
 ) throw()
{
  enzo_->SolveHydroEquations (data_block, enzo_->CycleNumber, dt);
}

