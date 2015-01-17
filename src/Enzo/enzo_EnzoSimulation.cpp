// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#include "simulation.hpp"


//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n)
{

  TRACE("EnzoSimulation::EnzoSimulation()");

  problem_ = new EnzoProblem;

  initialize();
}

//----------------------------------------------------------------------

EnzoSimulation::~EnzoSimulation() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulation::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  Simulation::pup(p);

  TRACEPUP;

  if (p.isUnpacking()) {
    EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			  field_descr());
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{

  // Initialize EnzoConfig parameters

  initialize_config_();

  // Call initialize() on base Simulation class
  Simulation::initialize();

  // Initialize EnzoBlock static variables
  EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			field_descr());
}

//----------------------------------------------------------------------

const Factory * EnzoSimulation::factory() const throw()
{ 
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize_config_() throw()
{
  if (config_ == NULL) config_ = new EnzoConfig;

  static_cast<EnzoConfig*>(config_)->read(parameters_);

}

//======================================================================

