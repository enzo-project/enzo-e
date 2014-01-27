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

extern CProxy_SimulationCharm proxy_simulation;

//----------------------------------------------------------------------

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

  monitor_ = Monitor::instance();

  monitor_->set_active (group_process->is_root());
  monitor_->header();
  monitor_->print ("","BEGIN ENZO-P");

  // Print initial baseline memory usage

  Memory * memory = Memory::instance();
  monitor_->print("Memory","bytes %lld bytes_high %lld",
		  memory->bytes(), memory->bytes_high());


  // open parameter file, displaying usage if invalid

  if (PARALLEL_ARGC != 2) {
    // Print usage if wrong number of arguments
    char buffer [ERROR_LENGTH];
    sprintf (buffer, "\nUsage: %s %s <parameter-file>\n\n", 
	     PARALLEL_RUN,PARALLEL_ARGV[0]);
    for (int i=0; i<PARALLEL_ARGC; i++) {
      PARALLEL_PRINTF("%d %s\n",i,PARALLEL_ARGV[i]);
    }
    ERROR("Main()",buffer);
  }

  const char * parameter_file = PARALLEL_ARGV[1];

  //--------------------------------------------------

  proxy_main     = thishandle;

  proxy_simulation = CProxy_EnzoSimulationCharm::ckNew
    (parameter_file, strlen(parameter_file)+1);

}

PARALLEL_MAIN_END


//======================================================================
#include "enzo.def.h"
//======================================================================
