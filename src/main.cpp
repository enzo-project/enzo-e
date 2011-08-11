// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file main.cpp

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "parallel.hpp"
#include "monitor.hpp"

#include "main.hpp"

//----------------------------------------------------------------------

CProxy_Main proxy_main;

void Main::p_exit(int count)
{
  TRACE("Main::p_exit");
  count_exit_++;
  if (count_exit_ >= count) {
    count_exit_ = 0;
    monitor_->print ("END ENZO-P");
    //    unit_finalize();
    // Fake unit_init() for index.php (test.hpp is not included since
    // enzo.ci and test.ci conflict)
    PARALLEL_PRINTF ("UNIT TEST END\n");
    PARALLEL_EXIT;
  }
}

#include "main.def.h"

#endif
