// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @brief    Implementation of main-level CHARM entry functions in main.ci

#include "cello.hpp"
#include "test.hpp"
#include "parallel.hpp"
#include "monitor.hpp"
#include "main.hpp"

//----------------------------------------------------------------------

CProxy_Main proxy_main;

#ifdef CHARM_ENZO
#include "simulation.hpp"
extern CProxy_Simulation proxy_simulation;
#include "enzo_finalize.hpp"
#endif

//----------------------------------------------------------------------

void Main::p_exit(int count)
{
  DEBUG("Main::p_exit");
  count_exit_++;
  unit_finalize();
  if (count_exit_ >= count) {
    count_exit_ = 0;
    exit_();
  }
}

void Main::exit_()
{

#ifdef CHARM_ENZO

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  if (simulation) {
    enzo_finalize(simulation);
  }

#endif

  if (Monitor::instance()) {
    Monitor::instance()->print ("","END CELLO");
  }
  PARALLEL_EXIT;
}


//----------------------------------------------------------------------

void Main::p_checkpoint(int count, std::string dir_name)
{
  
  count_checkpoint_++;
  if (count_checkpoint_ >= count) {
    count_checkpoint_ = 0;
    // Write parameter file

    char dir_char[255];
    strcpy(dir_char,dir_name.c_str());

    // --------------------------------------------------
    // ENTRY: #1 OutputCheckpoint::write_simulation()-> Simulation::s_write()
    // ENTRY: checkpoint if Simulation is root
    // --------------------------------------------------
#ifdef CHARM_ENZO
    CkCallback callback(CkIndex_Simulation::r_write_checkpoint(),proxy_simulation);
    CkStartCheckpoint (dir_char,callback);
#endif
  }
  // --------------------------------------------------
}


//----------------------------------------------------------------------

void Main::p_output_enter()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_output_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_output_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_output_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_enter()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_compute_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_continue()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_compute_continue();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_compute_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_enter()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_stopping_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_balance()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_stopping_balance();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_stopping_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_enzo_matvec()
{
#ifdef CHARM_ENZO
  CProxy_EnzoBlock * enzo_array = (CProxy_EnzoBlock*)proxy_simulation.ckLocalBranch()->hierarchy()->block_array();
  enzo_array->p_enzo_matvec();
#endif
}

//----------------------------------------------------------------------

void Main::p_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_enter()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_adapt_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_initial_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_initial_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_end()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_adapt_end();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_next()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_adapt_next();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_called()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_adapt_called();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_adapt_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_refresh_enter()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_refresh_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_refresh_exit()
{
#ifdef CHARM_ENZO
  proxy_simulation.ckLocalBranch()->hierarchy()->block_array()->p_refresh_exit();
#endif
}

//----------------------------------------------------------------------

#if defined(CHARM_ENZO)
#  include "main_enzo.def.h"
#elif defined(CHARM_SIMULATION)
#  include "main_simulation.def.h"
#elif defined(CHARM_MESH)
#  include "main_mesh.def.h"
#else
#  include "main.def.h"
#endif

