// See LICENSE_CELLO file for license and copyright information

/// @file     _performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    Private include file for the \ref Performance component

#ifndef _PERFORMANCE_HPP
#define _PERFORMANCE_HPP

//----------------------------------------------------------------------

/* Performance attributes */

enum attribute_enum {
  attribute_undefined, // 0 is an undefined attribute
  attribute_timestep,  // Simulation timesteps [monotonic]
  attribute_level,     // AMR hierarchy level
  attribute_component, // software component [memory]
  attribute_function,  // code function
  num_attributes = attribute_function // Number of attribute types
};

/* Performance counters */

enum counter_enum {
  counter_undefined, // 0 is an undefined counter
  counter_comm_send_bytes,     // Amount of data sent from this thread
  counter_comm_recv_bytes,     // Amount of data sent from this thread
  counter_comm_send_time,      // Time spent sending data
  counter_comm_recv_time,      // Time spent receiving data
  counter_comm_global_time,    // Time spent in collective communication
  counter_comm_send_count,     // Number of sends
  counter_comm_recv_count,     // Number of receives
  counter_comm_global_count,   // Number of barriers/reductions

#ifdef CONFIG_USE_PAPI
  counter_time_user,           // CPU time in user code of region
  counter_time_sys,            // CPU time in system of region
  counter_cpu_flop_count,      // Number of floating point operations
  mem_count,           // Number of memory accesses
#endif /* CONFIG_USE_PAPI */
  counter_time_real,           // Wallclock time of region
  counter_time_sim,            // Simulation time
  counter_mem_curr_bytes,      // Current number of bytes allocated
  counter_mem_high_bytes,      // Maximum number of bytes allocated
  counter_mem_new_count,       // Number of calls to allocate memory
  counter_mem_delete_count,    // Number of calls to deallocate memory
  counter_mem_new_bytes,       // Number of bytes allocated
  counter_mem_delete_bytes,    // Number of bytes deallocated
  counter_disk_read_bytes,     // Number of bytes read from disk
  counter_disk_write_bytes,    // Number of bytes written to disk
  counter_disk_read_time,      // Time spent reading from disk
  counter_disk_write_time,     // Time spent writing to disk
  counter_user_patch_count,    // Number of grid patches in each level
  counter_user_cell_count,     // Number of grid cells in each level
  counter_user_particle_count, // Number of particles
  num_counters = counter_user_particle_count
};

/* Performance functions */

enum function_enum {
  function_undefined, // 0 is an undefined function
  num_functions = function_undefined
};

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <vector>
#include <stack>
#include <string>
#include <sstream>
#include <sys/resource.h>

#ifdef __linux__
#   include <unistd.h>
#endif
#ifdef CONFIG_USE_PAPI
#  include "papi.h"
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "performance_Timer.hpp"
#include "performance_Counters.hpp"
#include "performance_Papi.hpp"
#include "performance_Performance.hpp"


#endif /* _PERFORMANCE_HPP */
