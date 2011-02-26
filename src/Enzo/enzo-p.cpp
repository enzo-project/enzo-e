// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @todo      support multiple input files
/// @brief     Cello main

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//**********************************************************************
#include PARALLEL_CHARM_INCLUDE(enzo_p.decl.h)
//**********************************************************************

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN

{

  // initialize

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  // create globals

  Monitor    * monitor    = new Monitor;
  Parameters * parameters = new Parameters (monitor);

  // only display output from root process

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

  EnzoSimulation * simulation = new EnzoSimulation (monitor);

  // read parameter file

  simulation->initialize(fp);

  // run the simulation

  simulation->run();

  // display footer text

  monitor->print ("END ENZO-P");

  // clean up

  delete simulation;
  delete parameters;
  delete monitor;

  // exit

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(enzo_p.def.h)
