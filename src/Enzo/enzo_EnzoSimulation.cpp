// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

/* #define CHECK_MEMORY */  /* calls mtrace() */

#define TRACE_PARAMETERS

#include "cello.hpp"

#include "enzo.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#ifdef CHECK_MEMORY
#   include <mcheck.h>
#endif

#include "simulation.hpp"

CProxy_EnzoSimulation proxy_enzo_simulation;

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char         parameter_file[],
 int                n,
 int                node_size) throw ()
  : CBase_EnzoSimulation(parameter_file, n),
    node_size_(node_size)
{
#ifdef CHECK_MEMORY
  mtrace();
#endif

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

  p | node_size_;

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

  const int ipn = CkMyPe();
  const int npn = MIN(node_size_,CkNumPes());

  if ((ipn % npn) == 0) {
#ifdef TRACE_PARAMETERS
    CkPrintf ("%d CALLING read_parameters(%d/%d)\n",CkMyPe(),ipn,npn);
    fflush(stdout);
#endif
    read_parameters_();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::read_parameters_()
{

  const int ipn = CkMyPe();
  const int npn = MIN(node_size_,CkNumPes());
#ifdef TRACE_PARAMETERS
  CkPrintf ("%d BEGIN read_parameters(%d/%d)\n",CkMyPe(),ipn,npn);
  fflush(stdout);
#endif
  //const int npn = CkNumPes();

  // Read parameter file
  parameters_ = new Parameters(parameter_file_.c_str(),monitor_);

#ifdef TRACE_PARAMETERS
  CkPrintf ("%d END read_parameters(%d/%d)\n",CkMyPe(),ipn,npn);
  fflush(stdout);
#endif
  // Then tell next Simulation object in node to read parameter file
  if (((ipn + 1) % npn) != 0 && (ipn < CkNumPes())) {
#ifdef TRACE_PARAMETERS
    CkPrintf ("%d CALLING p_read_parameters(%d/%d)\n",CkMyPe(),CkMyPe()+1,npn);
    fflush(stdout);
#endif
    proxy_enzo_simulation[CkMyPe()+1].p_read_parameters();
  }

  // Everybody synchronizes afterwards with barrier

#ifdef TRACE_PARAMETERS
  CkPrintf ("%d CALLING r_startup_finished()\n",CkMyPe());
  fflush(stdout);
#endif
  CkCallback callback (CkIndex_EnzoSimulation::r_startup_finished(NULL),
			 thisProxy);
  contribute(callback);

}

//----------------------------------------------------------------------

void EnzoSimulation::r_startup_finished (CkReductionMsg *msg)
{
#ifdef TRACE_PARAMETERS
  CkPrintf ("%d BEGIN r_startup_finished()\n",CkMyPe());
  fflush(stdout);
#endif
  delete msg;

  problem_ = new EnzoProblem;

  initialize_config_();

  initialize();
#ifdef TRACE_PARAMETERS
  CkPrintf ("%d END   r_startup_finished()\n",CkMyPe());
  fflush(stdout);
#endif
}

//----------------------------------------------------------------------

void EnzoSimulation::r_write_checkpoint()
{
  CkPrintf("EnzoSimulation::r_write_checkpoint() BEGIN\n");
  problem()->output_wait(this);
  CkPrintf("EnzoSimulation::r_write_checkpoint() END\n");
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

