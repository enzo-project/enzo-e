// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add ParallelCharm (for CkExit(), etc.)
/// @date     2009-10-16
/// @brief    [\ref Parallel] Interface for the Parallel class

#ifndef PARALLEL_PARALLEL_HPP
#define PARALLEL_PARALLEL_HPP

class Parallel {

  /// @class    Parallel
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Container group for hierarchical Parallel objects,
  /// e.g. GroupProcessMpi, GroupThreadOmp, etc.

public: // interface

  /// Initialize a Parallel object (singleton design pattern)
  Parallel(GroupProcess * process_group = 0,
	   GroupThread  * thread_group = 0) throw()
    : process_group_(process_group),
      thread_group_(thread_group)
  { };

  /// Get total number of processors
  GroupProcess * process_group() const throw()
  { return process_group_; };

  /// Get rank of this process
  GroupThread * thread_group() const throw()
  { return thread_group_; };

  /// Return whether this is the root or not
  bool is_root() const throw() 
  { return process_group_->is_root() && thread_group_->is_root(); };

private: // attributes

  /// Pointer to a group of distributed-memory processes
  GroupProcess * process_group_;

  /// Pointer to a group of shared-memory threads
  GroupThread * thread_group_;

};

#endif /* PARALLEL_PARALLEL_HPP */

