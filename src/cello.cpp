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

    //    printf ("\n");
    //    printf ("     The Laboratory for Computational Astrophysics proudly presents:\n");
    printf ("\n");
    printf ("    =================================================================\n");
    printf ("\n");
    printf ("    oooooooooooo                                          ooooo ooooo \n");
    printf ("    `888'     `8                                          `888' `888' \n");
    printf ("     888         ooo. .oo.     oooooooo  .ooooo.           888   888  \n");
    printf ("     888oooo8    `888P\"Y88b   d'\"\"7d8P  d88' `88b          888   888  \n");
    printf ("     888    \"     888   888     .d8P'   888   888 8888888  888   888  \n");
    printf ("     888       o  888   888   .d8P'  .P 888   888          888   888  \n");
    printf ("    o888ooooood8 o888o o888o d8888888P  `Y8bod8P'         o888o o888o\n");
    printf ("\n");
    printf ("    =================================================================\n");
    printf ("              E N Z O : T H E   N E X T  G E N E R A T I O N\n");
    printf ("    =================================================================\n");
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
