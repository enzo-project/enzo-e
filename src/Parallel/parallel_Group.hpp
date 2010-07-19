// $Id: parallel_group.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_HPP
#define PARALLEL_GROUP_HPP

/// @file     parallel_Group.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the ParallelGroup class

class ParallelGroup {

  /// @class    ParallelGroup
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for hierarchical parallelism

public: // interface

  /// Initialize the ParallelGroup object
  ParallelGroup(int process_min, int process_max,
		int thread_min,  int thread_max)
  {
    process_range_[0] = process_min;
    process_range_[1] = process_max;
    thread_range_[0] = thread_min;
    thread_range_[1] = thread_max;
  }

  /// Equality operator

  bool operator == (const ParallelGroup & group) throw()
  {
    return (process_range_[0] == group.process_range_[0] &&
	    process_range_[1] == group.process_range_[1] &&
	    thread_range_[0] == group.thread_range_[0] &&
	    thread_range_[1] == group.thread_range_[1]);

  }

  /// Return process rank

  void process_range (int * min, int * max) 
  { 
    *min = process_range_[0];
    *max = process_range_[1];
  }

  /// Return thread rank

  void thread_range (int * min, int * max) 
  { 
    *min = thread_range_[0];
    *max = thread_range_[1];
  };

private: // attributes

  int process_range_[2];
  int thread_range_[2];

};

#endif /* PARALLEL_GROUP_HPP */

