// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoSimulationMpi user-dependent class member functions

#ifndef CONFIG_USE_CHARM

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSimulationMpi::EnzoSimulationMpi
(
 const char * parameter_file,
 const GroupProcess * group_process ) throw ()
  : SimulationMpi(parameter_file,group_process)
{
  problem_ = new EnzoProblem;
}

//----------------------------------------------------------------------

EnzoSimulationMpi::~EnzoSimulationMpi() throw()
{
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::initialize() throw()
{
  initialize_config_();

  SimulationMpi::initialize();
  EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),field_descr());
}

//----------------------------------------------------------------------

const Factory * EnzoSimulationMpi::factory() const throw()
{ 
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//----------------------------------------------------------------------

void EnzoSimulationMpi::initialize_config_() throw()
{
  TRACE("BEGIN EnzoSimulationMpi::initialize_config_");
  if (config_ == NULL) config_ = new EnzoConfig;

  static_cast<EnzoConfig*>(config_)->read(parameters_);
  TRACE("END   EnzoSimulationMpi::initialize_config_");
}

#endif /* ! CONFIG_USE_CHARM */
