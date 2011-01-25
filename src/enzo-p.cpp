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

void usage(int argc, char ** argv);
void exit(Monitor *,GroupProcess *);

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN

{

  //==================================================
  // INITIALIZE
  //==================================================

  //-------------------------
  // initialize parallelism
  //-------------------------

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  //-------------------------
  // create globals
  //-------------------------

  Error   * error   = new Error;
  Monitor * monitor = new Monitor;

  //-------------------------
  // initialize monitor
  //-------------------------

  monitor->set_active(parallel->rank()==0);

    
  monitor->print ("BEGIN ENZO-P");

  monitor->header();

  //-------------------------
  // input parameters
  //-------------------------

  Parameters * parameters = new Parameters (monitor);

  FILE *fp = 0;

  // open parameter file

  if (PARALLEL_ARGC == 2) {
    fp = fopen(PARALLEL_ARGV[1],"r");
  }

  // print usage if invalid

  if ((PARALLEL_ARGC == 2 && !fp) || 
      (PARALLEL_ARGC != 2)) {
    if (parallel->is_root()) usage(PARALLEL_ARGC,PARALLEL_ARGV);
    PARALLEL_EXIT;
  }

  ASSERT ("enzo-p", "File pointer NULL", fp != 0);


  //-------------------------
  // Initialize simulation
  //-------------------------

  EnzoSimulation simulation (error,monitor);

  // read parameter file

  simulation.initialize(fp);
  
  //==================================================
  // FINALIZE
  //==================================================

  monitor->print ("END ENZO-P");

  delete parameters;
  delete monitor;
  delete error;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(enzo_p.def.h)

void usage(int argc, char ** argv)
{
  fprintf (stderr,"\nUsage: %s %s <parameter-file>\n\n",
	   PARALLEL_RUN,argv[0]);
}


