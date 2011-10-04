// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_COUNTER_HPP
#define PARALLEL_COUNTER_HPP

/// @file     parallel_Counter.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-04
/// @brief    [\ref Parallel] Declaration of the Counter class
///

class Counter {

  /// @class    Counter
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Class for terminating CHARM++ parallel "loops"

public: // interface

  private:
    int count_max_;
    int count_curr_;
  public:
    Counter (int count_max)
      : count_max_(count_max),
	count_curr_(0)
    {}
    int next ()
    {
      count_curr_++;
      if (count_curr_ >= count_max_) {
	count_curr_ = 0;
	
      return count_ = (count_max_ + count_ - 1) % count_max_;
    }
  };
}
#endif /* PARALLEL_COUNTER_HPP */

