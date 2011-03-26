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

int num_simulations;

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

    // display header text

    monitor->header();

    monitor->print ("BEGIN ENZO-P");

    // open parameter file, calling usage() if invalid

    if (PARALLEL_ARGC < 2) {
      // Print usage if wrong number of arguments
      char buffer [ERROR_LENGTH];
      sprintf (buffer,
	       "\nUsage: %s %s <parameter-file>\n\n", 
	       PARALLEL_RUN,PARALLEL_ARGV[0]);
      ERROR("main",buffer);
    }

    // Read in parameters

#ifdef CONFIG_USE_CHARM
    count_ = 0;
#endif

    int index_simulation;

    num_simulations = PARALLEL_ARGC - 1;

    for (index_simulation = 0; 
	 index_simulation < num_simulations; 
	 index_simulation++) {
      
      char * parameter_file = PARALLEL_ARGV[index_simulation + 1];

#ifdef CONFIG_USE_CHARM

      // If using CHARM, save the Main proxy, and create the
      // EnzoSimulationCharm group (one copy per process)

      mainProxy = thishandle;

      CProxy_EnzoSimulationCharm::ckNew(parameter_file, 
					strlen(parameter_file),
					index_simulation);

#else

      Simulation * simulation = 
	new EnzoSimulationMpi (parameter_file,group_process,index_simulation);

      ASSERT ("main()","Failed to create Simulation object",simulation != 0);

      // Initialize the simulation

      simulation->initialize();

      // Run the simulation

      simulation->run();

      delete simulation;
      
#endif

    } // for (int index_simulation ...)

#ifndef CONFIG_USE_CHARM    
    // display footer text

    Monitor::instance()->print ("END ENZO-P");

    // clean up

    delete group_process;

    // finalize unit testing

    unit_finalize();

    // exit

    PARALLEL_EXIT;
#endif

  };

#ifdef CONFIG_USE_CHARM
void enzo_exit(int index_simulation)
  {
    PARALLEL_PRINTF ("Simulation %d procs %d count %d\n",
		     index_simulation,
		     CkNumPes(),
		     count_);
    count_++;
    if (count_ == num_simulations * CkNumPes()) {
      Monitor::instance()->print ("END ENZO-P");
      unit_finalize();
      PARALLEL_EXIT;
    }
  };

private:
  int count_;
#endif

PARALLEL_MAIN_END


//======================================================================
#include PARALLEL_CHARM_INCLUDE(enzo.def.h)
//======================================================================
