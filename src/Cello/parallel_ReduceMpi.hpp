// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Reduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Mar 11 12:07:45 PST 2011
/// @brief    [\ref Parallel] Declaration and implementation of ReduceMpi

#ifndef PARALLEL_REDUCE_MPI_HPP
#define PARALLEL_REDUCE_MPI_HPP

#include "error.hpp"

class ReduceMpi : public Reduce {

  /// @class    ReduceMpi
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Implementation of ReduceMpi

public: // interface

  /// Constructor
  ReduceMpi(const GroupProcess * group_process) throw()
    : Reduce (group_process)
  { /* EMPTY */ };
    
  /// Destructor
  virtual ~ReduceMpi() throw()
  { /* EMPTY */ };

  /// Local reduction of the given integer value
  virtual int reduce_int
  (  int              local,
     enum_reduce_op   reduce_op )  throw();

  /// Local reduction of the given double value
  virtual double reduce_double
  (  double              local,
     enum_reduce_op   reduce_op )  throw();

};

#endif /* PARALLEL_REDUCE_MPI_HPP */

