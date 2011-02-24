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
void header(Monitor *);

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

  Monitor * monitor = new Monitor;

  //-------------------------
  // initialize monitor
  //-------------------------

  monitor->set_active(parallel->rank()==0);

    
  monitor->print ("BEGIN ENZO-P");

  header(monitor);

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

  EnzoSimulation simulation (monitor);

  // read parameter file

  simulation.initialize(fp);

  simulation.run();

  //==================================================
  // FINALIZE
  //==================================================

  monitor->print ("END ENZO-P");

  delete parameters;
  delete monitor;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(enzo_p.def.h)

//----------------------------------------------------------------------

void usage(int argc, char ** argv)
{
  fprintf (stderr,"\nUsage: %s %s <parameter-file>\n\n",
	   PARALLEL_RUN,argv[0]);
}


//----------------------------------------------------------------------

void header(Monitor * monitor)
{
  monitor->print ("");
  monitor->print ("    =================================================================");
  monitor->print ("");
  monitor->print ("    oooooooooooo                                          ooooooooo.   ");
  monitor->print ("    `888'     `8                                          `888   `Y88. ");
  monitor->print ("     888         ooo. .oo.     oooooooo  .ooooo.           888   .d88' ");
  monitor->print ("     888oooo8    `888P\"Y88b   d'\"\"7d8P  d88' `88b          888ooo88P'  ");
  monitor->print ("     888    \"     888   888     .d8P'   888   888 8888888  888         ");
  monitor->print ("     888       o  888   888   .d8P'  .P 888   888          888         ");
  monitor->print ("    o888ooooood8 o888o o888o d8888888P  `Y8bod8P'         o888o        ");
  monitor->print ("");
  monitor->print ("    =================================================================");
  monitor->print ("              E N Z O : T H E   N E X T  G E N E R A T I O N");
  monitor->print ("    =================================================================");
}
