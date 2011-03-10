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

#ifdef CONFIG_USE_MPI
  GroupProcess * parallel = GroupProcessMpi::create();
#else
  GroupProcess * parallel = GroupProcess::create();
#endif

  // initialize unit testing

  unit_init(parallel->rank(), parallel->size());

  // only display output from root process

  Monitor * monitor = Monitor::instance();

  monitor->set_active(parallel->rank()==0);

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
    if (parallel->is_root()) {
      fprintf (stderr,"\nUsage: %s %s <parameter-file>\n\n", 
	       PARALLEL_RUN,PARALLEL_ARGV[0]);
    }
    PARALLEL_EXIT;
  }

  // Read in parameters

  Parameters * parameters = new Parameters;

  parameters->read(fp); // MEMORY LEAK

  //--------------------------------------------------
  parameters->set_current_group ("Parallel");
  //--------------------------------------------------

  // parameter: Parallel::method

  std::string parallel_method = 
    parameters->list_value_string(0,"method","serial");

  Simulation * simulation = 0;
  if (parallel_method == "charm") {
    INCOMPLETE("main","EnzoSimulationCharm unfinished");
    //    simulation = new EnzoSimulationCharm (parameters);
  } else if (parallel_method == "serial") {
    simulation = new EnzoSimulationSerial (parameters);
  }

  ASSERT ("main()","Illegal Parallel::method parameter @@@\n",simulation != 0);

  // read parameter file

  simulation->initialize();

  // run the simulation

  simulation->run();

  // display footer text

  monitor->print ("END ENZO-P");

  // clean up

  delete simulation;
  delete parallel;
  delete parameters;

  // finalize unit testing

  unit_finalize();

  // exit

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(enzo_p.def.h)
