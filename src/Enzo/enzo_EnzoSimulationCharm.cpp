// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#include "simulation.hpp"


//----------------------------------------------------------------------

EnzoSimulationCharm::EnzoSimulationCharm
(
 const char         parameter_file[],
 int                n) throw ()
  : SimulationCharm(parameter_file, n)
{

  problem_ = new EnzoProblem;

// #ifdef CONFIG_USE_PROJECTIONS
//   traceRegisterUserEvent("Compute",10);
// #endif

  initialize();

}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::initialize() throw()
{
  SimulationCharm::initialize();
  EnzoBlock::initialize(parameters_,field_descr());
}

//----------------------------------------------------------------------

const Factory * EnzoSimulationCharm::factory() const throw()
{ 
  DEBUG("EnzoSimulationCharm::factory()");
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::run() throw()
{
  c_initial();
}

//======================================================================

#endif /* CONFIG_USE_CHARM */
