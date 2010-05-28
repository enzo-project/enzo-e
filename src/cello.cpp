// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      cello.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

//----------------------------------------------------------------------

#include <mpi.h>

#include "cello.hpp"

#include "error.hpp"
#include "parallel.hpp"
#include "schedule.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

//----------------------------------------------------------------------

void usage(int argc, char ** argv);

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{

  try {

    // INITIALIZE PARALLEL

    Parallel * parallel = Parallel::instance();

    parallel->initialize(&argc, &argv);

    // INITALIZE MONITOR

    Monitor * monitor = Monitor::instance();

    monitor->print ("CELLO BEGIN");

    monitor->header();

    // INPUT PARAMETERS

    FILE *fp = 0;
    if (argc == 2) {
      fp = fopen(argv[1],"r");
      if ( !fp ) usage(argc,argv);
    } else {
      usage(argc,argv);
    }

    ASSERT ("cello", "File pointer NULL", fp != 0);

    // READ PARAMETERS

    Parameters * parameters = Parameters::instance();

    parameters->read(fp);

    // INITIALIZE SIMULATION

    Schedule schedule;

    schedule.initialize_simulation();

    // RUN SIMULATION

    schedule.execute_simulation();

    // FINALIZE SIMULATION

    schedule.terminate_simulation();

    monitor->print ("CELLO END");

    parallel->finalize();

  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }

}

void usage(int argc, char ** argv)
{
  Parallel * parallel = Parallel::instance();

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
