// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Counters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-19 12:55:54
/// @brief    [\ref Performance] Interface for the Counters class

#ifndef PERFORMANCE_COUNTERS_HPP
#define PERFORMANCE_COUNTERS_HPP

class Counters {

  /// @class    Counters
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Represent counter values for a single
  /// set of attributes

public: // interface

  /// Initialize a Counters object
  Counters(size_t num_attributes, size_t num_counters)
    : num_attributes_(num_attributes),
      num_counters_(num_counters)
    {
      attributes_     = new int       [num_attributes];
      counters_start_ = new long long [num_counters];
      counters_stop_  = new long long [num_counters];
    }

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Counters()
    {
      delete [] attributes_;
      delete [] counters_start_;
      delete [] counters_stop_;
    }

  /// Copy constructor
  Counters(const Counters & classname) throw()
  {
  }

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p | num_attributes_;
    p | num_counters_;
    PUParray(p,attributes_,num_attributes_);
    PUParray(p,counters_start_,num_counters_);
    PUParray(p,counters_stop_,num_counters_);
    
  }
#endif

private: // attributes

  int num_attributes_;
  int num_counters_;

  /// Array of attribute values
  int       * attributes_;

  /// Array of the current counter values at the start of the region
  long long * counters_start_;

  /// Array of changes to the counter values
  long long * counters_stop_;

};

#endif /* PERFORMANCE_COUNTERS_HPP */

