// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      cello.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

//----------------------------------------------------------------------

#include <string>

#include "cello.hpp"

#include "parallel.hpp"
#include "simulation.hpp"

#include PARALLEL_CHARM_INCLUDE(cello.decl.h)

//----------------------------------------------------------------------

void usage(int argc, char ** argv);
void exit(Monitor *,GroupProcess *);

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN

{

  // Initialize parallelism

  PARALLEL_INIT;

  GroupProcess * parallel = GroupProcess::create();

  // INITALIZE "GLOBALS" (Parameters, Error, Monitor)

  Global * global = new Global;

  Monitor    * monitor    = global->monitor();
  Parameters * parameters = global->parameters();

  monitor->set_active(parallel->rank()==0);

    
  monitor->print ("CELLO BEGIN");

  monitor->header();

  // INPUT PARAMETERS

  FILE *fp = 0;
  if (PARALLEL_ARGC == 2) {
    fp = fopen(PARALLEL_ARGV[1],"r");
    if ( !fp ) {
      if (parallel->rank()==0) usage(PARALLEL_ARGC,PARALLEL_ARGV);
      PARALLEL_EXIT;
    }
  } else {
    if (parallel->rank()==0) usage(PARALLEL_ARGC,PARALLEL_ARGV);
    PARALLEL_EXIT;
  }

  ASSERT ("cello", "File pointer NULL", fp != 0);

  // READ PARAMETERS

  parameters->read(fp);

  // INITIALIZE SIMULATION

  Simulation simulation (global);

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END;

#include PARALLEL_CHARM_INCLUDE(cello.def.h)

void usage(int argc, char ** argv)
{
#ifdef CONFIG_USE_MPI
    fprintf (stderr,"\nUsage: mpirun [ options ] %s <parameter-file>\n\n",argv[0]);
#else
    fprintf (stderr,"\nUsage: %s <parameter-file>\n\n",argv[0]);
#endif
}


