// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_PARALLEL_GROUP_HPP
#define PARALLEL_PARALLEL_GROUP_HPP

/// @file     parallel_ParallelGroup.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:27 PDT 2010
/// @brief    Declaration of the ParallelGroup class

class ParallelGroup {

  /// @class    ParallelGroup
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for hierarchical parallelism

public: // interface

  /// Initialize the ParallelGroup object
  ParallelGroup(int size=1, int rank=0) throw() :
    size_(size),
    rank_(rank)
  { 
  }

  /// Number of compute elements in the ParallelGroup
  int size() throw()
  { return size_; };

  /// Rank of the compute element in the ParallelGroup
  int rank() throw()
  {  return rank_; };

  /// True iff rank() is 0
  bool is_root() throw()
  {  return rank_==0; };

  /// Synchronize between all compute elements in the ParallelGroup
  virtual void barrier() throw() { };

  /// Synchronize between two compute elements in the ParallelGroup
  virtual void sync(int rank) throw() { };

protected: // attributes

  /// Number of compute elements in the ParallelGroup
  int size_;

  /// Rank of this compute element in the ParallelGroup
  int rank_;

};

#endif /* PARALLEL_PARALLEL_GROUP_HPP */

