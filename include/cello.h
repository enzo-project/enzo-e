#ifndef CELLO_DEF
#define CELLO_DEF

/**********************************************************************
 * EDIT CONFIGURATION SETTINGS BELOW
 **********************************************************************/

#define CONFIG_PRECISION_SINGLE
/* #define CONFIG_PRECISION_DOUBLE */

#define CONFIG_USE_MPI
/* #define CONFIG_USE_PAPI */

#define USE_MEMORY

/*********************************************************************
 * PRECISION DECLARATIONS
 *********************************************************************/

#ifdef CONFIG_PRECISION_SINGLE
#   define Scalar float
#   define SCALAR_SCANF  "%f"
#   define SCALAR_PRINTF "%e "
#   define SCALAR_MPI     MPI_FLOAT
#   define SCALAR_HDF5    H5T_NATIVE_FLOAT
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#   define Scalar double
#   define SCALAR_SCANF  "%lf"
#   define SCALAR_PRINTF "%le "
#   define SCALAR_MPI     MPI_DOUBLE
#   define SCALAR_HDF5    H5T_NATIVE_DOUBLE
#endif

/*********************************************************************
 * PERFORMANCE DECLARATIONS
 **********************************************************************/

/* Performance attributes */

enum type_perf_attribute {
  attribute_undefined, // 0 is an undefined attribute
  attribute_timestep,  // Simulation timesteps [monotonic]
  attribute_level,     // AMR hierarchy level
  attribute_component, // software component [memory]
  attribute_function,  // code function
  num_attributes = attribute_function
};

/* Performance counters */

enum type_perf_counter {
  counter_undefined, // 0 is an undefined counter
#ifdef CONFIG_USE_MPI
  comm_send_bytes,     // Amount of data sent from this thread
  comm_recv_bytes,     // Amount of data sent from this thread
  comm_send_time,      // Time spent sending data
  comm_recv_time,      // Time spent receiving data
  comm_global_time,    // Time spent in collective communication
  comm_send_count,     // Number of sends
  comm_recv_count,     // Number of receives
  comm_global_count,   // Number of barriers/reductions
#endif /* CONFIG_USE_MPI */

#ifdef CONFIG_USE_PAPI
  time_user,           // CPU time in user code of region
  time_sys,            // CPU time in system of region
  cpu_flop_count,      // Number of floating point operations
  mem_count,           // Number of memory accesses
#endif /* CONFIG_USE_PAPI */
  time_real,           // Wallclock time of region
  time_sim,            // Simulation time
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
  num_counters = user_particle_count
};

/* Performance functions */

enum type_perf_functions {
  function_undefined, // 0 is an undefined function
  num_functions = function_undefined
};

/*********************************************************************
 * GLOBAL FUNCTIONS
 **********************************************************************/

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/*********************************************************************
 * COMPONENTS
 **********************************************************************/

enum type_component {
  component_undefined,
  component_amr,
  component_analysis,
  component_array,
  component_control,
  component_disk,
  component_error,
  component_field,
  component_memory,
  component_method,
  component_monitor,
  component_parallel,
  component_parameters,
  component_particles,
  component_performance,
  component_portal,
  component_problem,
  component_simulation,
  num_components
};


/* System includes */

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include <map>

#endif /* CELLO_DEF */
