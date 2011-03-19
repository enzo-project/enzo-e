// $Id: enzo_EnzoSimulationCharm.cpp 2116 2011-03-17 22:11:13Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Move timestep reductions into Timestep object
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 Parameters * parameters,
 GroupProcess * group_process) throw ()
  : EnzoSimulation(parameters,group_process)
{
}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

