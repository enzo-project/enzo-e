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

#define TRACE_RESTART


#ifdef TRACE_RESTART
#   define TRACE_RESTART_BLOCK(MSG,BLOCK)       \
  CkPrintf ("%d TRACE_RESTART %s %s \n",        \
            CkMyPe(),BLOCK->name().c_str(),     \
            std::string(MSG).c_str());          \
  fflush(stdout);
#else
#   define TRACE_RESTART_BLOCK(MSG,BLOCK) /* ... */
#endif

// SEE enzo_control_restart.cpp
