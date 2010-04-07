// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      cello.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

#include "cello.h"

#include "parallel.hpp"
#include "schedule.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

void usage(int argc, char ** argv) 
{
#ifdef CONFIG_USE_MPI
  fprintf (stderr,"Usage: mpirun [ options ] %s <parameter-file>\n\n",argv[0]);
#else
  fprintf (stderr,"Usage: %s <parameter-file>\n\n",argv[0]);
#endif
  exit(1);
}

int main(int argc, char ** argv)
{

  try {

    printf ("\n");   
    printf ("        Astrophysics and Cosmology for Supercomputers:\n");
    printf ("\n");   
    printf ("           .oooooo.             oooo  oooo            \n");
    printf ("          d8P'  `Y8b            `888  `888            \n");
    printf ("         888           .ooooo.   888   888   .ooooo.  \n");
    printf ("         888          d88' `88b  888   888  d88' `88b \n");
    printf ("         888          888ooo888  888   888  888   888 \n");
    printf ("         `88b    ooo  888    .o  888   888  888   888 \n");
    printf ("          `Y8bood8P'  `Y8bod8P' o888o o888o `Y8bod8P' \n");
    printf ("\n");
    printf ("        ==============================================\n");   
    printf ("        E N Z O : T H E   N E X T  G E N E R A T I O N\n");
    printf ("        ==============================================\n");   
    printf ("\n");

    // INITIALIZE PARALLEL

    ParallelMpi mpi;

    mpi.initialize(&argc, &argv);

    // INPUT PARAMETERS

    FILE *fp = 0;
    if (argc == 2) {
      fp = fopen(argv[1],"r");
      if ( !fp ) usage(argc,argv);
    } else {
      usage(argc,argv);
    }

    assert (fp != 0);

    // INITALIZE MONITOR

    Monitor monitor;
    monitor.print ("CELLO BEGIN");

    // INITIALIZE SCHEDULE

    Schedule schedule(&monitor);

    schedule.read_parameters(fp);

    // INITIALIZE SIMULATION

    schedule.initialize_simulation();

    // RUN SIMULATION

    schedule.execute_simulation();

    // FINALIZE SIMULATION

    schedule.terminate_simulation();

    monitor.print ("CELLO END");

    mpi.finalize();

  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }

}
