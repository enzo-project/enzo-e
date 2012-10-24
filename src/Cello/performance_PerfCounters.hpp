// See LICENSE_CELLO file for license and copyright information

/// @file     performance_PerfCounters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-19 12:55:54
/// @brief    [\ref Performance] Interface for the PerfCounters class

#ifndef PERFORMANCE_COUNTERS_HPP
#define PERFORMANCE_COUNTERS_HPP

class PerfCounters {

  /// @class    PerfCounters
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Represent counter values for a single
  /// region

public: // interface

  /// Initialize a PerfCounters object
  PerfCounters(size_t num_counters)
    : num_counters_(num_counters)
    {
      counters_start_ = new long long [num_counters];
      counters_stop_  = new long long [num_counters];
    }

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~PerfCounters()
    {
      delete [] counters_start_;
      delete [] counters_stop_;
    }

  /// Copy constructor
  PerfCounters(const PerfCounters & classname) throw()
  {
  }

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | num_counters_;
    PUParray(p,counters_start_,num_counters_);
    PUParray(p,counters_stop_,num_counters_);
    
  }
#endif

private: // attributes

  int num_counters_;

  /// Array of the current counter values at the start of the region
  long long * counters_start_;

  /// Array of changes to the counter values
  long long * counters_stop_;

};

#endif /* PERFORMANCE_COUNTERS_HPP */

