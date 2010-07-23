// $Id: parallel.hpp 1275 2010-03-09 01:05:14Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Parallel.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Add ParallelCharm (for CkExit(), etc.)
/// @date     2009-10-16
/// @brief    Interface for the Parallel class

#ifndef PARALLEL_PARALLEL_HPP
#define PARALLEL_PARALLEL_HPP

class Parallel {

  /// @class    Parallel
  /// @ingroup  Parallel
  /// @brief Container group for hierarchical Parallel objects,
  /// e.g. GroupProcessMpi, GroupThreadOmp, etc.

public: // interface

  /// Initialize a Parallel object (singleton design pattern)
  Parallel(GroupProcess * process_group = 0,
	   GroupThread  * thread_group = 0)
    : process_group_(process_group),
      thread_group_(thread_group)
  { };

  /// Initialize
  void initialize(int * argc = 0, char ***argv = 0)
  {
    if (process_group_) process_group_->initialize(argc,argv);
    if (thread_group_)  thread_group_->initialize();
  }

  /// Finalize
  void finalize()
  {
    if (process_group_) process_group_->finalize();
    if (thread_group_)  thread_group_->finalize();
  }

  /// Get total number of processors
  GroupProcess * process_group() { return process_group_; };

  /// Get rank of this process
  GroupThread * thread_group() { return thread_group_; };

private: // attributes

  /// Pointer to a group of distributed-memory processes
  GroupProcess * process_group_;

  /// Pointer to a group of shared-memory threads
  GroupThread * thread_group_;

};

#endif /* PARALLEL_PARALLEL_HPP */

