// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of SimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#include "simulation.hpp"

//----------------------------------------------------------------------

SimulationCharm::SimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n),
    block_sync_(0)
{
  TRACE("SimulationCharm::SimulationCharm");

  // derived class should call initialize()
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
  TRACE("SimulationCharm::~SimulationCharm()");
}

//----------------------------------------------------------------------

// void SimulationCharm::run() throw()
// {
//   TRACE("SimulationCharm::run()");
//   initial();
// }

//----------------------------------------------------------------------
//======================================================================

#endif /* CONFIG_USE_CHARM */
