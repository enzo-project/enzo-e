// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Reduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Mar 11 12:07:39 PST 2011
/// @brief    [\ref Parallel] Declaration and implementation of ReduceSerial

#ifndef PARALLEL_REDUCE_SERIAL_HPP
#define PARALLEL_REDUCE_SERIAL_HPP

#include "error.hpp"

class ReduceSerial : public Reduce {

  /// @class    Reduce
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Implementation of ReduceSerial

public: // interface

  /// Constructor
  ReduceSerial(GroupProcess * group_process) throw()
    : Reduce (group_process)
  { /* EMPTY */ };
    
  /// Destructor
  ~ReduceSerial() throw()
  { /* EMPTY */ };

  /// Parallel reduction of an integer
  virtual int reduce_int
  (  int            value,
     enum_reduce_op reduce_op )  throw()
  { return value; }

  /// Parallel reduction of a double
  virtual double reduce_double
  (  double         value,
     enum_reduce_op reduce_op )  throw()
  { return value; }

};

#endif /* PARALLEL_REDUCE_SERIAL_HPP */

