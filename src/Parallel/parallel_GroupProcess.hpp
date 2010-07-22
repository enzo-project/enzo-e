// $Id: parallel_GroupProcess.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_GROUP_PROCESS_HPP
#define PARALLEL_GROUP_PROCESS_HPP

/// @file     parallel_GroupProcess.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 16:40:35 PDT 2010
/// @brief    Declaration of the GroupProcess class

class GroupProcess : public Group {

  /// @class    GroupProcess
  /// @ingroup  Parallel
  /// @brief    Group of distributed memory processes

public: // interface

  /// Initialize the GroupProcess object
  GroupProcess(int process_first, int process_count)
    : Group(process_count),
      process_first_(process_first)
  {
  }

private: // attributes

  int process_first_;

};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

