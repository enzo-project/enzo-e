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
  Affinity(GroupProcess * processes = 0,
	   GroupThread  * threads   = 0) throw()
    : processes_(processes),
      threads_(threads),
      process_rank_(processes ? processes->rank() : 0),
      process_size_(processes ? processes->size() : 1)
  { }

  /// Equality operator

  bool operator == (const Affinity & affinity) throw()
  {
    return (process_rank() == affinity.process_rank() &&
	    thread_rank()  == affinity.thread_rank());
  }

  /// Return process rank

  int process_rank () throw()
  { return process_rank_; };

  /// Return process rank

  int process_size () throw()
  { return process_size_; };

  /// Return thread rank

  int thread_rank () throw()
  { return threads_ ? threads_->rank() : 0; };

  /// Return thread rank

  int thread_size () throw()
  { return threads_ ? threads_->size() : 0; };

  /// Return the Parallel class for distributed parallelism

  const GroupProcess * processes() throw()
  { return processes_; };

  /// Return the Parallel class for threaded parallelism

  const GroupThread * threads() throw()
  { return threads_; };

  /// Whether this is the root process / thread

  bool is_root() throw()
  {return process_rank() == 0 && thread_rank() == 0; } ;

private: // attributes

  GroupProcess * processes_;
  GroupThread * threads_;

  int process_rank_;
  int process_size_;

};

#endif /* PARALLEL_AFFINITY_HPP */

