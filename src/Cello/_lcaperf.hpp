// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-14
/// @brief    Include file for the \ref Lcaperf component
///
///   LcaPerf()
///   ~LcaPerf()
///   initialize()	 initialize
///   begin()		 set_active
///   end()		 
///                       register
///                       update
///
/// Attributes
///
/// [create attribute set]
///
///   new_attribute ()	 attribute_new
///   attribute()	 attribute_set
///   delete_attribute()	 attribute_del
///
/// Counters
///
/// [create counter set]
///
///   new_counter()	 counter_new
///   increment()	 counter_inc
///   assign()		 counter_set
///   delete_counter()	 counter_del
///   value()		 counter_val
///   print()             counters_print
///
/// Regions
///
/// [activate/deactivate counter set for regions]
/// [associate attribute sets to regions]
///
///   start()		 region_start
///   stop()      	 region_stop

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include <map>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <limits>
#include <sstream>

#ifdef CONFIG_USE_MPI
#   include <mpi.h>
#endif
#ifdef CONFIG_USE_PAPI
#  include <papi.h>
#endif

#include "lcaperf.def"
#include "lcaperf_attributes.hpp"
#include "lcaperf_counters.hpp"
#include "lcaperf_counters_basic.hpp"
#include "lcaperf_counters_user.hpp"
#include "lcaperf_counters_deriv.hpp"
#ifdef CONFIG_USE_PAPI
#include "lcaperf_counters_papi.hpp"
#endif
#ifdef CONFIG_USE_MPI
#include "lcaperf_counters_mpi.hpp"
#endif
#ifdef CONFIG_TRACE_MEM
#include "lcaperf_counters_mem.hpp"
#endif
#include "lcaperf_it_counter_keys.hpp"

#include "lcaperf_lcaperf.hpp"
