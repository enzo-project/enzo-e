// $Id: parallel_Group.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_HPP
#define PARALLEL_GROUP_HPP

/// @file     parallel_Group.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the Group class

class Group {

  /// @class    Group
  /// @ingroup  Parallel
  /// @brief    Generalization of MPI rank for hierarchical parallelism

public: // interface

  /// Initialize the Group object
  Group(int size=1, int rank=0):
    size_(size),
    rank_(rank)
  { };

  /// Number of compute elements in the Group
  int size() { return size_; };

  /// Rank of the compute element in the Group
  int rank() { return rank_; };

  /// Synchronize between all compute elements in the Group
  virtual void barrier() { };

  /// Synchronize between two compute elements in the Group
  virtual void wait() { };

private: // attributes

  /// Number of compute elements in the Group
  int size_;

  /// Rank of the compute element in the Group
  int rank_;

};

#endif /* PARALLEL_GROUP_HPP */

