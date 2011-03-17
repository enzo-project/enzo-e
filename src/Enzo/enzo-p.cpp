// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @todo      support multiple input files
/// @brief     Cello main

//----------------------------------------------------------------------

#include "test.hpp"

#include "enzo.hpp"

//**********************************************************************
#include PARALLEL_CHARM_INCLUDE(enzo_p.decl.h)
//**********************************************************************

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN

{

  // initialize parallel

  PARALLEL_INIT;

  // Create global parallel process group object
  GroupProcess * group_process = GroupProcess::create();

  // initialize unit testing

  int rank = group_process->rank();
  int size = group_process->size();

  unit_init(rank, size);

  Monitor * monitor = Monitor::instance();

  // only display output from root process
  monitor->set_active(rank == 0);

  // display header text

  monitor->header();

  monitor->print ("BEGIN ENZO-P");

  // open parameter file, calling usage() if invalid

  FILE *fp = 0;

  if (PARALLEL_ARGC == 2) {
    fp = fopen(PARALLEL_ARGV[1],"r");
  }
  if ((PARALLEL_ARGC == 2 && !fp) || 
      (PARALLEL_ARGC != 2)) {
    if (group_process->is_root()) {
      fprintf (stderr,"\nUsage: %s %s <parameter-file>\n\n", 
	       PARALLEL_RUN,PARALLEL_ARGV[0]);
    }
    PARALLEL_EXIT;
  }

  // Read in parameters

  Parameters * parameters = new Parameters;

  parameters->read(fp);

  Simulation * simulation = 0;

#ifdef CONFIG_USE_CHARM
  simulation = new EnzoSimulationCharm (parameters, group_process);
#else
  simulation = new EnzoSimulationMpi (parameters,group_process);
#endif

  ASSERT ("main()","Failed to create Simulation object",simulation != 0);

  // Initialize the simulation

  simulation->initialize();

  // Run the simulation

  simulation->run();

  // display footer text

  monitor->print ("END ENZO-P");

  // clean up

  delete simulation;
  delete group_process;
  delete parameters;

  // finalize unit testing

  unit_finalize();

  // exit

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(enzo_p.def.h)
