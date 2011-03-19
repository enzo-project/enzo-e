// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationMpi user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationMpi::EnzoSimulationMpi
(
 Parameters * parameters,
 GroupProcess * group_process) throw ()
  : EnzoSimulation(parameters,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationMpi::~EnzoSimulationMpi() throw()
{
}
