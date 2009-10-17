#ifndef CELLO_DEF
#define CELLO_DEF

/*********************************************************************
* Define either CONFIG_PRECISION_SINGLE or CONFIG_PRECISION_DOUBLE
**********************************************************************/

#define CONFIG_PRECISION_SINGLE
/* #define CONFIG_PRECISION_DOUBLE */

/*********************************************************************
* GLOBAL DECLARATIONS
**********************************************************************/

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
#   define SCALAR_MPI    MPI_DOUBLE
#   define SCALAR_HDF5    H5T_NATIVE_DOUBLE
#endif

/* Performance attribute values for Performance::attr_component */

enum type_perf_group {
  perf_group_amr,
  perf_group_analysis,
  perf_group_array,
  perf_group_control,
  perf_group_disk,
  perf_group_error,
  perf_group_field,
  perf_group_memory,
  perf_group_method,
  perf_group_monitor,
  perf_group_parallel,
  perf_group_parameters,
  perf_group_particles,
  perf_group_performance,
  perf_group_portal,
  perf_group_problem,
  perf_group_simulation
};

/* System includes */

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include <map>

/* Include "global" components */

#include "error.hpp"
#include "performance.hpp"


#endif /* CELLO_DEF */
