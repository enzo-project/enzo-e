// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

/* #define CHECK_MEMORY */  /* calls mtrace() */

// #define TRACE_PARAMETERS
// #define DEBUG_ENZO_SIMULATION

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
 int                n)
  : BASE_ENZO_SIMULATION(parameter_file, n)
{
#ifdef CHECK_MEMORY
  mtrace();
#endif

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION EnzoSimulation()\n",CkMyPe());
  fflush(stdout);
#endif  

  // Synchronize to ensure all EnzoSimulation objects exist before
  // reading parameters

  CkCallback callback (CkIndex_EnzoSimulation::r_startup_begun(NULL),
                       thisProxy);
  contribute(callback);

}

//----------------------------------------------------------------------

EnzoSimulation::~EnzoSimulation()
{
}

//----------------------------------------------------------------------

void EnzoSimulation::pup (PUP::er &p)
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION EnzoSimulation::pup()\n",CkMyPe());
  fflush(stdout);
#endif  
  // NOTE: change this function whenever attributes change

  BASE_ENZO_SIMULATION::pup(p);

  TRACEPUP;

  if (p.isUnpacking()) {
    EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			  field_descr());
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::r_startup_begun (CkReductionMsg *msg)
{

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif  

  delete msg;

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif  

  problem_ = new EnzoProblem;

  initialize_config_();

  initialize();
#ifdef TRACE_PARAMETERS
  CkPrintf ("%d END   r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif
}

//----------------------------------------------------------------------

void EnzoSimulation::r_write_checkpoint()
{
  performance_->start_region(perf_output);
  monitor()->print ("Output","begin writing checkpoint...");
  problem()->output_wait(this);
  monitor()->print ("Output","end writing checkpoint");
  performance_->stop_region(perf_output);
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize_config_() throw()
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION begin initialize_config()\n",CkMyPe());
  fflush(stdout);
#endif  

  config_ = static_cast<EnzoConfig*>(&g_enzo_config);

  // char buffer[40];
  // sprintf (buffer,"config-%02d.text",CkMyPe());
  // FILE * fp = fopen (buffer,"w");
  // EnzoConfig * enzo_config = static_cast<EnzoConfig*>(config_);
  //  enzo_config->write(fp);
  // fclose(fp);
  
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION end initialize_config()\n",CkMyPe());
  fflush(stdout);
#endif  

}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION initialize()\n",CkMyPe());
  fflush(stdout);
#endif  
  // Call initialize() on base Simulation class
  Simulation::initialize();

  // Initialize EnzoBlock static variables
  EnzoBlock::initialize(static_cast<EnzoConfig*>(config_),
			field_descr());

}

//----------------------------------------------------------------------

const Factory * EnzoSimulation::factory() const throw()
{ 
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION factory()\n",CkMyPe());
  fflush(stdout);
#endif  
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//======================================================================

