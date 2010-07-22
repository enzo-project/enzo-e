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
    : Group(process_max_plus - process_min)
  {
    process_range_[0] = process_min;
    process_range_[1] = process_max;
  }

  /// Return process rank

  void process_range (int * min, int * max) 
  { 
    *min = process_range_[0];
    *max = process_range_[1];
  }

  /// Return thread rank

private: // attributes

  int process_range_[2];

};

#endif /* PARALLEL_GROUP_PROCESS_HPP */

