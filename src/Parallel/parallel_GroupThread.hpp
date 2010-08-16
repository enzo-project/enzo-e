// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_THREAD_HPP
#define PARALLEL_GROUP_THREAD_HPP

/// @file     parallel_GroupThread.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul 22 12:36:17 PDT 2010
/// @brief    Declaration of the GroupThread class

class GroupThread : public ParallelGroup {

  /// @class    GroupThread
  /// @ingroup  Parallel
  /// @brief    ParallelGroup of shared memory threads

public: // interface

  /// Initialize the GroupThread object
  GroupThread(int size = 1, int rank = 0)
    : ParallelGroup(size,rank)
  {  }

private: // attributes


};

#endif /* PARALLEL_GROUP_THREAD_HPP */

