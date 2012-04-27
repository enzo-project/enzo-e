#ifndef LCAPERF_H
#define LCAPERF_H

//======================================================================
//
// File:        lcaperf.hpp
//
// Description: Parallel performance monitoring class
//
//----------------------------------------------------------------------
//
// Classes:     LcaPerf
//
//======================================================================

//   LcaPerf()
//   ~LcaPerf()
//   initialize()	 initialize
//   begin()		 set_active
//   end()		 
//                       register
//                       update
//
// Attributes
//
// [create attribute set]
//
//   new_attribute ()	 attribute_new
//   attribute()	 attribute_set
//   delete_attribute()	 attribute_del
//
// Counters
//
// [create counter set]
//
//   new_counter()	 counter_new
//   increment()	 counter_inc
//   assign()		 counter_set
//   delete_counter()	 counter_del
//   value()		 counter_val
//   print()             counters_print
//
// Regions
//
// [activate/deactivate counter set for regions]
// [associate attribute sets to regions]
//
//   start()		 region_start
//   stop()      	 region_stop

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

#ifdef CONFIG_USE_MPI
#   include <mpi.h>
#endif
#ifdef CONFIG_USE_PAPI
#  include <papi.h>
#endif

namespace lca {

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

//======================================================================
// LcaPerf base class
//======================================================================
class LcaPerf {


 public: // functions

  /// LcaPerf constructor
  LcaPerf();

  /// LcaPerf destructor
  ~LcaPerf();

  /// Add a region name: used for printing
  void new_region (const char * region_name);

  /// Create a new attribute of the given type
  void new_attribute (const char * attribute_name,
		      int          attribute_type);

  /// Return the value of the attribute of the given type
  void attribute (const char        * attribute_name,
		  const void        * attribute_value,
		  int                 attribute_type);

  /// Delete the named attribute
  void delete_attribute (const char * attribute_name);

  /// Create a new counter of the given type
  void new_counter  (const char *      counter_name,
		     counter_type counter_type);

  /// Increment the given user counter
  void increment (const char * counter_name,
		  long long    value);

  /// Assign the value to the user counter
  void assign    (const char * counter_name,
		  long long    value);

  /// Delete the named counter of the given type
  void delete_counter (const char * counter_name);

  /// Initialize lcaperf
  void initialize (const char * filename = 0);

  /// Finalize lcaperf
  void finalize ();

  /// Begin collecting performance data in files with the give suffix
  void begin ();

  /// Stop collecting performance data in files with the given suffix
  void end   ();

  /// Return the counter value for the currently assigned attribute
  long long value (const char * set,
		   const char * key,
		   const char * counter);

  /// Identify the beginning of a code region
  void start     (const char * region_name);

  /// Identify the end of a code region
  void stop      (const char * region_name);

  //--------------------------------------------------

  /// Print the header for subsequent print()'s
  void header ();

  /// Print all keys and associated counter values
  void print ();

  //----------------------------------------------------------------------

 protected: // attributes

  /// Attributes object
  Attributes     attributes_;

  /// List of Counter object pointers
  std::map<std::string,Counters *> counters_;

  /// List of Regions encountered
  std::vector<std::string> regions_;

};
}
#endif
