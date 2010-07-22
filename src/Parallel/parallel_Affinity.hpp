// $Id: parallel_Affinity.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_AFFINITY_HPP
#define PARALLEL_AFFINITY_HPP

/// @file     parallel_Affinity.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the Affinity class

class Affinity {

  /// @class    Affinity
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for hierarchical parallelism

public: // interface

  /// Initialize the Affinity object
  Affinity(int process_rank = 0,
	   int thread_rank  = 0) throw()
    : process_rank_ (process_rank),
      thread_rank_ (thread_rank),
      processes_(0),
      threads_(0),
      group_(0)
  {}

  /// Equality operator

  bool operator == (const Affinity & affinity) throw()
  {
    return (process_rank_ == affinity.process_rank_ &&
	    thread_rank_  == affinity.thread_rank_);
  }

  /// Return process rank

  int process_rank () { return process_rank_; };

  /// Return thread rank

  int thread_rank () { return process_rank_; };

  /// Return the Parallel class for distributed parallelism

  const GroupProcess * processes() { return processes_; };

  /// Return the Parallel class for threaded parallelism

  const GroupThread * threads() { return threads_; };

  /// Return the Parallel group

  const ParallelGroup * group() { return group_; };

private: // attributes

  int process_rank_;
  int thread_rank_;

  GroupProcess * processes_;
  GroupThread * threads_;

};

#endif /* PARALLEL_AFFINITY_HPP */

