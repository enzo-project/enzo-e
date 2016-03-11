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

CProxy_EnzoSimulation proxy_enzo_simulation;

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char         parameter_file[],
 int                n) throw ()
  : Simulation(parameter_file, n)
{
  // Synchronize to ensure all EnzoSimulation objects exist before
  // reading parameters

  CkCallback callback (CkIndex_EnzoSimulation::r_startup_begun(NULL),
			 thisProxy);
  contribute(callback);

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

void EnzoSimulation::r_startup_begun (CkReductionMsg *msg)
{
  delete msg;

  // Serialize reading parameters within each logical node
  if (CkMyRank() == 0) {
    read_parameters_();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::read_parameters_()
{
  const int ipn = CkMyRank();
  const int npn = CkNodeSize(ipn);

  // Read parameter file
  parameters_ = new Parameters(parameter_file_.c_str(),monitor_);

  // Then tell next Simulation object in node to read parameter file
  if (((ipn + 1) % npn) != 0) {
    proxy_enzo_simulation[CkMyPe()+1].p_read_parameters();
  }

  // Everybody synchronizes afterwards with barrier

  CkCallback callback (CkIndex_EnzoSimulation::r_startup_finished(NULL),
			 thisProxy);
  contribute(callback);
}

//----------------------------------------------------------------------

void EnzoSimulation::r_startup_finished (CkReductionMsg *msg)
{
  delete msg;

  problem_ = new EnzoProblem;

  initialize_config_();

  initialize();
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize_config_() throw()
{
  if (config_ == NULL) {
    config_ = new EnzoConfig;
  }

  static_cast<EnzoConfig*>(config_)->read(parameters_);

}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{

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

//======================================================================

