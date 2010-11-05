// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_PARALLEL_AFFINITY_HPP
#define PARALLEL_PARALLEL_AFFINITY_HPP

/// @file     parallel_ParallelAffinity.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the ParallelAffinity class

class ParallelAffinity {

  /// @class    ParallelAffinity
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for hierarchical parallelism

public: // interface

  /// Initialize the ParallelAffinity object
  ParallelAffinity(GroupProcess * processes = 0,
	   GroupThread  * threads   = 0) throw()
    : process_rank_(processes ? processes->rank() : 0),
      process_size_(processes ? processes->size() : 1),
      thread_rank_(threads ? threads->rank() : 0),
      thread_size_(threads ? threads->size() : 0)
  { }

  /// Equality operator

  bool operator == (const ParallelAffinity & affinity) const throw()
  {
    return (process_rank() == affinity.process_rank() &&
	    thread_rank()  == affinity.thread_rank());
  }

  /// Whether this is the root process / thread

  bool is_root() const throw()
  {return process_rank() == 0 && thread_rank() == 0; } ;

  /// Return process rank

  int process_rank () const throw()
  { return process_rank_; };

  /// Return process rank

  int process_size () const throw()
  { return process_size_; };

  /// Return thread rank

  int thread_rank () const throw()
  { return thread_rank_; };

  /// Return thread rank

  int thread_size () const throw()
  { return thread_size_; };

private: // attributes

  /// Rank in the process group
  int process_rank_;

  /// Size of the process group
  int process_size_;

  /// Rank in the thread group
  int thread_rank_;
  /// Size of the thread group
  int thread_size_;

};

#endif /* PARALLEL_PARALLEL_AFFINITY_HPP */

