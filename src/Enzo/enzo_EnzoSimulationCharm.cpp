// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulationCharm user-dependent class member functions

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

  TRACE("EnzoSimulationCharm::EnzoSimulationCharm()");

  problem_ = new EnzoProblem;

  initialize();
}

//----------------------------------------------------------------------

EnzoSimulationCharm::~EnzoSimulationCharm() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  SimulationCharm::pup(p);
  if (p.isUnpacking()) {
    EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			  field_descr());
  }
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::initialize() throw()
{
  
  initialize_config_();

  SimulationCharm::initialize();
  EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			field_descr());
}

//----------------------------------------------------------------------

const Factory * EnzoSimulationCharm::factory() const throw()
{ 
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//----------------------------------------------------------------------

void EnzoSimulationCharm::initialize_config_() throw()
{
  if (config_ == NULL) config_ = new EnzoConfig;

  static_cast<EnzoConfig*>(config_)->read(parameters_);

}

//======================================================================

