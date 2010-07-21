// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PERFORMANCE_COUNTERS_HPP
#define PERFORMANCE_COUNTERS_HPP

/// @file     performance_Counters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-19 12:55:54
/// @brief    Interface for the Counters class

class Counters {

  /// @class    Counters
  /// @ingroup  Performance
  /// @brief    Represent counter values for a single set of attributes

public: // interface

  /// Initialize a Counters object
  Counters(size_t num_attributes, size_t num_counters)
    {
      //      Memory * memory = Memory::instance();
      //      memory->begin_group(component_performance);
      a_  = new int       [num_attributes];
      c_  = new long long [num_counters];
      dc_ = new long long [num_counters];
      //      memory->end_group(component_performance);
    }

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Counters()
    {
      //      Memory * memory = Memory::instance();
      //      memory->begin_group(component_performance);
      delete [] a_;
      delete [] c_;
      delete [] dc_;
      //      memory->end_group(component_performance);
    }

  /// Copy constructor
  Counters(const Counters & classname) throw()
  {
    INCOMPLETE_MESSAGE("Counters::Counters","");
  }

  /// Assignment operator
  Counters & operator= (const Counters & classname) throw()
  {
    INCOMPLETE_MESSAGE("Counters::operator =","");
    return *this;
  }

private: // attributes

  /// Array of attribute values
  int       * a_;

  /// Array of the current counter values at the start of the region
  long long * c_;

  /// Array of changes to the counter values
  long long * dc_;

};

#endif /* PERFORMANCE_COUNTERS_HPP */

