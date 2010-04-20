// $Id: parallel_affinity.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef AFFINITY_HPP
#define AFFINITY_HPP

/// @file     parallel_affinity.hpp
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
      thread_rank_ (thread_rank)
  {}

  /// Equality operator

  bool operator == (const Affinity & affinity) throw()
  {
    return (process_rank_ == affinity.process_rank_ &&
	    thread_rank_  == affinity.thread_rank_);
  }

private: // attributes

  int process_rank_;
  int thread_rank_;

};

#endif /* AFFINITY_HPP */

