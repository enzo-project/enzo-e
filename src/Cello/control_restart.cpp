// See LICENSE_CELLO file for license and copyright information

/// @file     control_restart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions.
/// @ingroup  Control
///
/// This file controls adaptive mesh refinement on a distributed
/// array of octrees.

//--------------------------------------------------

#include "simulation.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"
#include "main.hpp"

#define TRACE_RESTART


#ifdef TRACE_RESTART
#   define TRACE_RESTART_BLOCK(MSG,BLOCK)       \
  CkPrintf ("TRACE_RESTART %s %s \n",           \
            BLOCK->name().c_str(),              \
            std::string(MSG).c_str());          \
  fflush(stdout);
#   define TRACE_RESTART_MAIN(MSG)              \
  CkPrintf ("TRACE_RESTART %s \n",              \
            std::string(MSG).c_str());          \
  fflush(stdout);
#else
#   define TRACE_RESTART_BLOCK(MSG,BLOCK) /* ... */
#   define TRACE_RESTART_MAIN(MSG) /* ... */
#endif

//----------------------------------------------------------------------

void Block::restart_enter_()
{
  TRACE_RESTART_BLOCK("restart_enter_",this);
  const std::string restart_dir  = cello::config()->initial_restart_dir;
  const std::string restart_file = cello::config()->initial_restart_file;
  if (index_.is_root()) {
    proxy_main.p_restart_enter(restart_dir,restart_file);
  }
}

//----------------------------------------------------------------------

void Main::p_restart_enter
(std::string restart_dir, std::string restart_file)
{
  // open hierarchy file
  TRACE_RESTART_MAIN("Main::restart_enter_ open hierarchy file");
  CkPrintf ("DEBUG_RESTART restart_dir = %s\n",restart_dir.c_str());
  CkPrintf ("DEBUG_RESTART restart_file = %s\n",restart_file.c_str());
  // read hierarchy file
  TRACE_RESTART_MAIN("Main::restart_enter_ read hierarchy file");
  // create block_array
  TRACE_RESTART_MAIN("Main::restart_enter_ create block array");
  // create IoEnzoReader array
  TRACE_RESTART_MAIN("Main::restart_enter_ create IoEnzoReader array");
  // initialize sync_file(num_io_reader)
  //  for (i_f = files in restart) {
  //    io_reader[i_f].insert(file_block);
  //  }
  // close hierarcy file
  TRACE_RESTART_MAIN("restart_enter_ close hierarchy file");
}
