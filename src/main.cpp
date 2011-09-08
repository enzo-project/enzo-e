// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @brief    Implementation of main-level CHARM entry functions in main.ci
/// @todo     Add callback function for cleanup, including UNIT TEST END, etc.

#ifdef CONFIG_USE_CHARM

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
#endif

void Main::p_exit(int count)
{
  count_exit_++;
  if (count_exit_ >= count) {
    count_exit_ = 0;
#ifdef CHARM_ENZO

    // DUPLICATE CODE IN src/enzo-p.cpp !!!

    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Parameters * parameters = simulation->parameters();

    int    cycle_final = parameters->value_integer("Testing:cycle_final",0);

    if (cycle_final != 0) {
      unit_assert (simulation->cycle()==cycle_final);
      monitor_->print ("cycle:  %d",simulation->cycle());
    }

    double time_final  = parameters->value_float("Testing:time_final",0.0);

    if (time_final != 0.0) {
      unit_assert (simulation->time()==time_final);
      monitor_->print ("time:  %g",simulation->time());
    }

#endif

    monitor_->print ("END CELLO");
    //    unit_finalize();
    // Fake unit_init() for index.php (test.hpp is not included since
    // enzo.ci and test.ci conflict)
    PARALLEL_PRINTF ("UNIT TEST END\n");
    PARALLEL_EXIT;
  }
}


#if defined(CHARM_ENZO)
#  include "main_enzo.def.h"
#elif defined(CHARM_SIMULATION)
#  include "main_simulation.def.h"
#elif defined(CHARM_MESH)
#  include "main_mesh.def.h"
#else
#  include "main.def.h"
#endif

#endif
