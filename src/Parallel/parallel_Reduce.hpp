// $Id: parallel_Reduce.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_REDUCE_HPP
#define PARALLEL_REDUCE_HPP

/// @file     parallel_Reduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Mar 11 12:07:34 PST 2011
/// @brief    [\ref Parallel] Interface for the Reduce class

#include "error.hpp"

class Reduce {

  /// @class    Reduce
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Implementation of reduction operations

public: // interface

  /// Constructor
  Reduce(GroupProcess * group_process) throw()
    : group_process_(group_process)
  { /* EMPTY */ };
    
  /// Destructor
  ~Reduce() throw()
  { /* EMPTY */ };

  /// Parallel int reduction of the stored local value
  virtual int reduce_int 
  (  int              value,
     enum_reduce_op   reduce_op ) throw()= 0;

  /// Parallel double reduction of the stored local value
  virtual double reduce_double 
  (  double           value,
     enum_reduce_op   reduce_op) throw()= 0;

protected: // attributes

  GroupProcess * group_process_;

};

#endif /* PARALLEL_REDUCE_HPP */

