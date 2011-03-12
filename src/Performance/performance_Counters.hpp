// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PERFORMANCE_COUNTERS_HPP
#define PERFORMANCE_COUNTERS_HPP

/// @file     performance_Counters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-19 12:55:54
/// @brief    [\ref Performance] Interface for the Counters class

class Counters {

  /// @class    Counters
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Represent counter values for a single
  /// set of attributes

public: // interface

  /// Initialize a Counters object
  Counters(size_t num_attributes, size_t num_counters)
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
    INCOMPLETE("Counters::Counters");
  }

  /// Assignment operator
  Counters & operator= (const Counters & classname) throw()
  {
    INCOMPLETE("Counters::operator =");
    return *this;
  }

private: // attributes

  /// Array of attribute values
  int       * attributes_;

  /// Array of the current counter values at the start of the region
  long long * counters_start_;

  /// Array of changes to the counter values
  long long * counters_stop_;

};

#endif /* PERFORMANCE_COUNTERS_HPP */

