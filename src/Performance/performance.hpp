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
    time_real_rel,       // Wallclock time of region
    first_counter = time_real_rel,
    time_real_abs,       // Elapsed wallclock time
    time_user_rel,       // CPU time in user code of region
    time_user_abs,       // Elapsed CPU time in user code
    time_sys_rel,        // CPU time in system of region
    time_sys_abs,        // Elapsed CPU time in system
    time_sim_abs,        // Simulation time
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
    cpu_mem_count,       // Number of memory accesses
    last_counter = cpu_mem_count
  };

  enum type_attr {
    attr_timestep,  // Simulation timesteps [monotonic]
    first_attr = attr_timestep,
    attr_level,     // AMR hierarchy level
    attr_component, // software component [memory]
    attr_function,  // code function
    attr_process,   // MPI process id
    last_attr = attr_process
  };

  enum type_metric {
    time_balance, // Relative time load balance at given level
    first_metric = time_balance,
    mem_balance,  // Relative memory load balance at given level
    comm_balance, // Relative communication load balance at given level
    last_metric = comm_balance
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
