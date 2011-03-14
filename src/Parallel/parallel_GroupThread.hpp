// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_THREAD_HPP
#define PARALLEL_GROUP_THREAD_HPP

/// @file     parallel_GroupThread.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:17 PDT 2010
/// @brief    [\ref Parallel] Declaration of the GroupThread class

class GroupThread {

  /// @class    GroupThread
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Group of shared memory threads

public: // interface

  /// Initialize the GroupThread object
  GroupThread(int size = 1, int rank = 0)
    : size_(size),
      rank_(rank)
  {  }

  /// Number of compute elements in the GroupThread
  int size() throw()
  { return size_; };

  /// Rank of the compute element in the GroupThread
  int rank() throw()
  {  return rank_; };

  /// True iff rank() is 0
  bool is_root() throw()
  {  return rank_==0; };

  /// Synchronize between all compute elements in the GroupThread
  virtual void barrier() throw() { };

  /// Synchronize between two compute elements in the GroupThread
  virtual void sync(int rank) throw() { };

private: // attributes

  /// Number of compute elements in the GroupThread
  int size_;

  /// Rank of this compute element in the GroupThead
  int rank_;

};

#endif /* PARALLEL_GROUP_THREAD_HPP */

