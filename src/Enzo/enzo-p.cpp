// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @todo      support multiple input files
/// @brief     Cello main

//----------------------------------------------------------------------

#include "test.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
  {

    // initialize parallel

    PARALLEL_INIT;

    // Create global parallel process group object
    GroupProcess * group_process = GroupProcess::create();

    // initialize unit testing

    int rank = group_process->rank();
    int size = group_process->size();

    unit_init(rank, size);

    Monitor * monitor = Monitor::instance();

    // only display output from root process
    monitor->set_active(rank == 0);

    // display header text

    monitor->header();

    monitor->print ("BEGIN ENZO-P");

    // open parameter file, calling usage() if invalid

    if (PARALLEL_ARGC != 2) {
      // Print usage if wrong number of arguments
      char buffer [ERROR_LENGTH];
      sprintf (buffer,
	       "\nUsage: %s %s <parameter-file>\n\n", 
	       PARALLEL_RUN,PARALLEL_ARGV[0]);
      ERROR("main",buffer);
    }

    // Read in parameters

    const char * parameter_file = PARALLEL_ARGV[1];

#ifdef CONFIG_USE_CHARM

    // If using CHARM, save the Main proxy, and create the
    // EnzoSimulationCharm group (one copy per process)

    mainProxy = thishandle;
    count_ = 0;

    CProxy_EnzoSimulationCharm::ckNew(parameter_file, strlen(parameter_file));

#else

    Simulation * simulation = 
      new EnzoSimulationMpi (parameter_file,group_process);

    ASSERT ("main()","Failed to create Simulation object",simulation != 0);

    // Initialize the simulation

    simulation->initialize();

    // Run the simulation

    simulation->run();

    // display footer text

    Monitor::instance()->print ("END ENZO-P");

    // clean up

    delete simulation;
    delete group_process;

    // finalize unit testing

    unit_finalize();

    // exit

    PARALLEL_PRINTF ("main exiting\n");

    PARALLEL_EXIT;

#endif

  };

#ifdef CONFIG_USE_CHARM
  void enzo_exit(void)
  {
    PARALLEL_PRINTF ("count_=%d\n",count_);
    count_++;
    if (count_ == CkNumPes()) {
      PARALLEL_PRINTF ("main exiting\n");
      CkExit();
    }
  };

private:
  int count_;
#endif

PARALLEL_MAIN_END


//======================================================================
#include PARALLEL_CHARM_INCLUDE(enzo.def.h)
//======================================================================
