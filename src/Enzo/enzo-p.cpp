// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @brief     Cello main
///
/// \mainpage Enzo-P / Cello
///
/// <a href="http://cello-project.org/">Cello</a> is an
/// object-oriented adaptive mesh refinement (AMR) software framework
/// for high performance scientific applications.  The framework is
/// scalable, easy to use, and portable across systems ranging from
/// laptops and PC's to the largest HPC systems available, including
/// Blue Waters, the National Science Foundation's Cray petascale
/// supercomputer at the University of Illinois at Urbana-Champaign.
/// Cello's mesh refinement uses the highly-scalable
/// "forest-of-octrees" approach, and the Charm++ parallel programming
/// system enables its high parallel scalability.
///
/// Development of Cello is driven by Enzo, a parallel computational
/// astrophysics and cosmology application. The goal is to efficiently
/// map Enzo's multi-resolution multi-physics capabilities onto large
/// parallel computers with potentially millions of computational
/// units. This "petascale" incarnation of Enzo being built on the
/// Cello framework is called Enzo-P.

//----------------------------------------------------------------------

#define CHARM_ENZO

#include "test.hpp"
#include "enzo.hpp"
#include "main.hpp"

#include "charm_enzo.hpp"

//----------------------------------------------------------------------

extern CProxy_EnzoSimulation proxy_enzo_simulation;
extern CProxy_Simulation proxy_simulation;

//----------------------------------------------------------------------
PARALLEL_MAIN_BEGIN
{

  // Initialize parallelization

  PARALLEL_INIT;

#ifdef PNG_1_2_X
  CkPrintf ("PNG_1_2_X\n");
#endif
#ifdef PNG_1_3_X
  CkPrintf ("PNG_1_3_X\n");
#endif
#ifdef PNG_1_4_X
  CkPrintf ("PNG_1_4_X\n");
#endif
#ifdef PNG_1_5_X
  CkPrintf ("PNG_1_5_X\n");
#endif
  

  // Check parameter file

  if (PARALLEL_ARGC != 2) {
    // Print usage if wrong number of arguments
   printf ("\nUsage: %s %s <parameter-file> [ +balancer <load-balancer> ]\n\n", 
	     PARALLEL_RUN,PARALLEL_ARGV[0]);
    p_exit(1);
  }
  
  const char * parameter_file = PARALLEL_ARGV[1];

  g_parameters.read(PARALLEL_ARGV[1]);
  g_parameters.write("parameters.out");
  g_enzo_config.read(&g_parameters);

  // Initialize unit testing

  const int ip = CkMyPe();
  const int np = CkNumPes();

  unit_init(ip,np);

  // Initialize Monitor

  monitor_ = Monitor::instance();
  monitor_->set_mode (monitor_mode_root);
  monitor_->header();
  monitor_->print ("","BEGIN ENZO-P");

  // Print initial baseline memory usage

  Memory * memory = Memory::instance();
  monitor_->print("Memory","bytes %lld bytes_high %lld",
		  memory->bytes(), memory->bytes_high());

#ifdef CONFIG_USE_PAPI
  int retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval > 0) {
    WARNING("Papi::init","PAPI library version mismatch!");
  } else if (retval < 0) {
    WARNING("Papi::init","PAPI initialization error!");
  }
#endif

 //--------------------------------------------------

  proxy_main     = thishandle;

  // --------------------------------------------------
  // ENTRY: #1 Main::Main() -> EnzoSimulation::EnzoSimulation()
  // ENTRY: create
  // --------------------------------------------------
  proxy_simulation = proxy_enzo_simulation = CProxy_EnzoSimulation::ckNew
    (parameter_file, strlen(parameter_file)+1);
  // --------------------------------------------------

}

PARALLEL_MAIN_END


//======================================================================
#include "enzo.def.h"
//======================================================================
