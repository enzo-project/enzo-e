// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

//----------------------------------------------------------------------

#define CHARM_ENZO

#include "test.hpp"

#include "enzo.hpp"
#include "main.hpp"

#include "charm_enzo.hpp"

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
extern CProxy_SimulationCharm proxy_simulation;
#endif

#ifndef CONFIG_USE_CHARM
#   include "enzo_finalize.hpp"
#endif

PARALLEL_MAIN_BEGIN
{

  // initialize parallelization

  PARALLEL_INIT;

  // Create global parallel process group object
  const GroupProcess * group_process = GroupProcess::create();

  // initialize unit testing

  int ip = group_process->rank();
  int np = group_process->size();

  unit_init(ip,np);

#ifdef CONFIG_USE_CHARM
  monitor_ = Monitor::instance();
#else
  Monitor * monitor_ = Monitor::instance();
#endif

  monitor_->set_active (group_process->is_root());
  monitor_->header();
  monitor_->print ("","BEGIN ENZO-P");

  // Print initial baseline memory usage

  Memory * memory = Memory::instance();
  monitor_->print("Memory","bytes %lld bytes_high %lld",
		  memory->bytes(), memory->bytes_high());


  // open parameter file, calling usage() if invalid

  if (PARALLEL_ARGC != 2) {
    // Print usage if wrong number of arguments
    char buffer [ERROR_LENGTH];
    sprintf (buffer, "\nUsage: %s %s <parameter-file>\n\n", 
	     PARALLEL_RUN,PARALLEL_ARGV[0]);
    ERROR("main",buffer);
  }

  // Read in parameters

  char * parameter_file = PARALLEL_ARGV[1];

  //--------------------------------------------------

#ifdef CONFIG_USE_CHARM

  proxy_main     = thishandle;

  proxy_simulation = CProxy_EnzoSimulationCharm::ckNew
    (parameter_file, strlen(parameter_file)+1);
  TRACE1("proxy_simulation = %p",&proxy_simulation);

  //--------------------------------------------------

#else /* ! CONFIG_USE_CHARM */

  Simulation * simulation = 
    new EnzoSimulationMpi (parameter_file,group_process);

  simulation->initialize();

  simulation->run();

  enzo_finalize(simulation);

  delete simulation;
  delete group_process;

  unit_finalize();

  PARALLEL_EXIT;

#endif
  //--------------------------------------------------
}

PARALLEL_MAIN_END


//======================================================================
#ifdef CONFIG_USE_CHARM
#  include "enzo.def.h"
#endif
//======================================================================
