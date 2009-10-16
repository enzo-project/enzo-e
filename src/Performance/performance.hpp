#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

#include "cello.h"

#include "performance_timer.hpp"

/** 
 *********************************************************************
 *
 * @file      performance.hpp
 * @brief     Class for collecting and allowing access to performance data
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Wed Oct 14 23:40:13 PDT 2009
 * @bug       
 * @note      
 *
 * Class for collecting and allowing access to performance data
 *
 * $Id$
 *
 *********************************************************************
 */

class Performance {

  enum type_counter {
    first_counter,       // Wallclock time of region
    time_real = first_counter,
    time_user,           // CPU time in user code of region
    time_sys,            // CPU time in system of region
    time_sim,            // Simulation time
    comm_send_bytes,     // Amount of data sent from this thread
    comm_recv_bytes,     // Amount of data sent from this thread
    comm_send_time,      // Time spent sending data
    comm_recv_time,      // Time spent receiving data
    comm_global_time,    // Time spent in collective communication
    comm_send_count,     // Number of sends
    comm_recv_count,     // Number of receives
    comm_global_count,   // Number of barriers/reductions
    mem_curr_bytes,      // Current number of bytes allocated
    mem_high_bytes,      // Maximum number of bytes allocated
    mem_new_count,       // Number of calls to allocate memory
    mem_delete_count,    // Number of calls to deallocate memory
    mem_new_bytes,       // Number of bytes allocated
    mem_delete_bytes,    // Number of bytes deallocated
    disk_read_bytes,     // Number of bytes read from disk
    disk_write_bytes,    // Number of bytes written to disk
    disk_read_time,      // Time spent reading from disk
    disk_write_time,     // Time spent writing to disk
    user_patch_count,    // Number of grid patches in each level
    user_cell_count,     // Number of grid cells in each level
    user_particle_count, // Number of particles
    cpu_flop_count,      // Number of floating point operations
    last_counter,        // Number of memory accesses
    cpu_mem_count = last_counter
  };

  enum type_attr {
    first_attr,          // Simulation timesteps [monotonic]
    attr_timestep = first_attr,
    attr_level,          // AMR hierarchy level
    attr_component,      // software component [memory]
    attr_function,       // code function
    last_attr,           // MPI process id
    attr_process = last_attr
  };

  enum type_metric {
    first_metric,        // Relative time load balance at given level
    time_balance = first_metric,
    mem_balance,         // Relative memory load balance at given level
    last_metric,         // Relative communication load balance at given level
    comm_balance = last_metric
  };

/** 
 *********************************************************************
 *
 * @class     Performance
 * @brief     
 * @ingroup   
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// 
  Performance();

  /// 
  ~Performance();

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// 


};

#endif /* PERFORMANCE_HPP */
