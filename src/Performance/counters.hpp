// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef COUNTERS_HPP
#define COUNTERS_HPP

/// @file     counters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-19 12:55:54
/// @brief    Interface for the Counters class

#include "cello.h"
#include "memory.hpp"

class Counters {

  /// @class    Counters
  /// @ingroup  Performance
  /// @brief    Represent counter values for a single set of attributes

public: // interface

  /// Initialize a Counters object
  Counters(size_t num_attributes, size_t num_counters)
    {
      Memory::begin_group(component_performance);
      a_  = new int       [num_attributes];
      c_  = new long long [num_counters];
      dc_ = new long long [num_counters];
      Memory::end_group(component_performance);
    }

  /// Delete a Counters object
  ~Counters()
    {
      Memory::begin_group(component_performance);
      delete [] a_;
      delete [] c_;
      delete [] dc_;
      Memory::end_group(component_performance);
    }


private: // attributes

  /// Array of attribute values
  int       * a_;

  /// Array of the current counter values at the start of the region
  long long * c_;

  /// Array of changes to the counter values
  long long * dc_;

};

#endif /* COUNTERS_HPP */

