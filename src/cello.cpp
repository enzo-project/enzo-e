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

void usage(int argc, char ** argv, Parallel * parallel);

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{

  try {

    // INITIALIZE PARALLEL

    ParallelCreate parallel_create;
    Parallel * parallel = parallel_create.create(parallel_mpi);

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
      if ( !fp ) usage(argc,argv,parallel);
    } else {
      usage(argc,argv,parallel);
    }

    ASSERT ("cello", "File pointer NULL", fp != 0);

    // READ PARAMETERS

    parameters->read(fp);

    // INITIALIZE SIMULATION

    Simulation simulation (global);

    monitor->print ("CELLO END");

    parallel->finalize();

  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }

}

void usage(int argc, char ** argv, Parallel * parallel)
{
  if (parallel->is_root()) {
#ifdef CONFIG_USE_MPI
    fprintf (stderr,"Usage: mpirun [ options ] %s <parameter-file>\n\n",argv[0]);
#else
    fprintf (stderr,"Usage: %s <parameter-file>\n\n",argv[0]);
#endif
  }
  parallel->finalize();
  exit(1);
}
