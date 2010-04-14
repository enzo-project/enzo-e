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

int main(int argc, char ** argv)
{

  try {

    // INITIALIZE PARALLEL

    Parallel * parallel = Parallel::instance();

    parallel->initialize(&argc, &argv);

    // INITALIZE MONITOR

    Monitor * monitor = Monitor::instance();

    monitor->print ("CELLO BEGIN");

    //    monitor->print ("");
    //    monitor->print ("     The Laboratory for Computational Astrophysics proudly presents:");
    monitor->print ("");
    monitor->print ("    =================================================================");
    monitor->print ("");
    monitor->print ("    oooooooooooo                                          ooooo ooooo ");
    monitor->print ("    `888'     `8                                          `888' `888' ");
    monitor->print ("     888         ooo. .oo.     oooooooo  .ooooo.           888   888  ");
    monitor->print ("     888oooo8    `888P\"Y88b   d'\"\"7d8P  d88' `88b          888   888  ");
    monitor->print ("     888    \"     888   888     .d8P'   888   888 8888888  888   888  ");
    monitor->print ("     888       o  888   888   .d8P'  .P 888   888          888   888  ");
    monitor->print ("    o888ooooood8 o888o o888o d8888888P  `Y8bod8P'         o888o o888o");
    monitor->print ("");
    monitor->print ("    =================================================================");
    monitor->print ("              E N Z O : T H E   N E X T  G E N E R A T I O N");
    monitor->print ("    =================================================================");
    monitor->print ("");

    // INPUT PARAMETERS

    FILE *fp = 0;
    if (argc == 2) {
      fp = fopen(argv[1],"r");
      if ( !fp ) usage(argc,argv);
    } else {
      usage(argc,argv);
    }

    assert (fp != 0);

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
