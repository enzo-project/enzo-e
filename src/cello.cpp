// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      cello.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

//----------------------------------------------------------------------

#include <string>

#include <mpi.h>

#include "cello.hpp"

#include "parallel.hpp"
#include "simulation.hpp"
#include "global.hpp"

//----------------------------------------------------------------------

void usage(int argc, char ** argv);
void exit(Monitor *,Parallel *);

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{

  try {

    // INITIALIZE PARALLEL

    parallel->initialize(&argc, &argv);

    // INITALIZE "GLOBALS" (Parameters, Error, Monitor)

    Global * global = new Global(parallel);

    Monitor    * monitor    = global->monitor();
    Parameters * parameters = global->parameters();
    
    monitor->print ("CELLO BEGIN");

    monitor->header();

    // INPUT PARAMETERS

    FILE *fp = 0;
    if (argc == 2) {
      fp = fopen(argv[1],"r");
      if ( !fp ) {
	if (parallel->is_root()) usage(argc,argv);
	exit(monitor,parallel);
      }
    } else {
      if (parallel->is_root()) usage(argc,argv);
      exit(monitor,parallel);
    }

    ASSERT ("cello", "File pointer NULL", fp != 0);

    // READ PARAMETERS

    parameters->read(fp);

    // INITIALIZE SIMULATION

    Simulation simulation (global);

    exit(monitor,parallel);
  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
    exit(1);
  }
}

void usage(int argc, char ** argv)
{
#ifdef CONFIG_USE_MPI
    fprintf (stderr,"\nUsage: mpirun [ options ] %s <parameter-file>\n\n",argv[0]);
#else
    fprintf (stderr,"\nUsage: %s <parameter-file>\n\n",argv[0]);
#endif
}

void exit(Monitor * monitor, Parallel * parallel)
{
  monitor->print ("CELLO END");
  parallel->finalize();
  exit(0);
}

