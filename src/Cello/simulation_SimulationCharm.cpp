// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of SimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "simulation.hpp"

#include "simulation_charm.hpp"
#include "mesh_charm.hpp"

#include "simulation.hpp"


//----------------------------------------------------------------------

SimulationCharm::SimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n)
{
  // derived class should call initialize()
}

//----------------------------------------------------------------------

SimulationCharm::~SimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void SimulationCharm::initialize() throw()
{
  Simulation::initialize();
}

//----------------------------------------------------------------------

void SimulationCharm::run() throw()
{
  c_initial();
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
