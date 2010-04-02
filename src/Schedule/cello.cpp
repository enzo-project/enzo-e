/// @file      cello.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main

#include "cello.h"

#include "parallel.hpp"h
#include "schedule.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

void usage(int argc, char ** argv) 
{
  fprintf (stderr,"Usage: %s <filename>\n",argv[0]);
  exit(1);
}

int main(int argc, char ** argv)
{

  try {

    printf ("\n	 ====================================================\n");
    printf ("		      ENZO: THE NEXT GENERATION\n\n");
    printf ("	  A parallel astrophysics and cosmology application\n\n");
    printf ("	     See CELLO_LICENSE for full license agreement\n");
    printf ("	 ====================================================\n\n");

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
