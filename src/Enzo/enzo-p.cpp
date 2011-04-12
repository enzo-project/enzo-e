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

#define MAX_OUTPUT 10 /* Maximum number of outputs going on at a time */

#ifdef CONFIG_USE_CHARM
CProxy_Main                proxy_main;
CProxy_EnzoSimulationCharm proxy_simulation;
#endif

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

  TRACE("");
  monitor_ = new Monitor;
  TRACE("");

  monitor_->set_active (true);
  TRACE("");

  // display header text
  TRACE("");

  monitor_->header();
  TRACE("");

  monitor_->print ("BEGIN ENZO-P");

  TRACE("");
  // open parameter file, calling usage() if invalid

  if (PARALLEL_ARGC != 2) {
    // Print usage if wrong number of arguments
    char buffer [ERROR_LENGTH];
    sprintf (buffer,
	     "\nUsage: %s %s <parameter-file>\n\n", 
	     PARALLEL_RUN,PARALLEL_ARGV[0]);
    ERROR("main",buffer);
  }

  TRACE("");
  // Read in parameters

  //--------------------------------------------------
#ifdef CONFIG_USE_CHARM
  proxy_main     = thishandle;

  // Clear counts
  count_exit_        = 0;
  count_prepare_        = 0;
  for (int i=0; i<MAX_OUTPUT; i++) {
    count_output_open_[i]  = 0;
    count_output_close_[i] = 0;
  }
#endif
  //--------------------------------------------------
  TRACE("");

     
  char * parameter_file = PARALLEL_ARGV[1];

  //--------------------------------------------------

#ifdef CONFIG_USE_CHARM

  // If using CHARM, create the EnzoSimulationCharm groups

  TRACE("");
  proxy_simulation = CProxy_EnzoSimulationCharm::ckNew
    (parameter_file, strlen(parameter_file)+1, 0);

  TRACE("");
  //--------------------------------------------------

#else /* ! CONFIG_USE_CHARM */

  Simulation * simulation = 
    new EnzoSimulationMpi (parameter_file,group_process, 0);

  ASSERT ("main()","Failed to create Simulation object",simulation != 0);

  // Initialize the simulation

  simulation->initialize();

  // Run the simulation

  simulation->run();

  // Delete the simulation

  delete simulation;
      
#endif

  //--------------------------------------------------
#ifndef CONFIG_USE_CHARM    
  // display footer text

  monitor_->print ("END ENZO-P");

  // clean up

  delete group_process;

  // finalize unit testing

  unit_finalize();

  // exit

  PARALLEL_EXIT;
#endif
  //--------------------------------------------------

};

//--------------------------------------------------
#ifdef CONFIG_USE_CHARM

//----------------------------------------------------------------------

void p_exit(int count)
{
  count_exit_++;
  if (count_exit_ >= count) {
    count_exit_ = 0;
    monitor_->print ("END ENZO-P");
    unit_finalize();
    PARALLEL_EXIT;
  }
};

//----------------------------------------------------------------------

void p_prepare(int count)
{
  count_prepare_++;
  if (count_prepare_ >= count) {
    count_prepare_ = 0;
    TRACE("main::p_prepare()");
    unit_finalize();
    PARALLEL_EXIT;
  }
};

//----------------------------------------------------------------------

//  --- Open output file and and initialize output data ---

void p_output_open(int count, int index, int cycle, double time)
{
  count_output_open_[index]++;
  if (count_output_open_[index] == count) {
    count_output_open_[index] = 0;
    // Request blocks to contribute index'th output data
    proxy_simulation.p_output(index,cycle,time);
    
  }
};

//----------------------------------------------------------------------

//  --- Accumulate block output contributions and write output to disk ---

void p_output_close(int count)
{
};

//----------------------------------------------------------------------

private:

int count_exit_;
int count_prepare_;
int count_output_open_[MAX_OUTPUT];
int count_output_close_[MAX_OUTPUT];

Monitor * monitor_;

#endif
//--------------------------------------------------

PARALLEL_MAIN_END


//======================================================================
#include PARALLEL_CHARM_INCLUDE(enzo.def.h)
//======================================================================
