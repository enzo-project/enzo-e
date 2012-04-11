// See LICENSE_CELLO file for license and copyright information

/// @file     charm_Counter.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-12
/// @brief    [\ref Charm] Declaration of the CHARM Counter class

#ifndef CHARM_COUNTER_HPP
#define CHARM_COUNTER_HPP

#ifdef CONFIG_USE_CHARM
#include "charm++.h"
class Counter {

  /// @class    Counter
  /// @ingroup  Charm
  /// @brief    [\ref Charm] 

public:
  Counter (int count_max) throw()
    : count_max_(count_max),
      count_curr_(0)
  {}

  int remaining() throw ()
  {
    if (count_max_ > 0) {
      count_curr_ = (count_max_ + count_curr_ - 1) % count_max_;  
    }
    return count_curr_;
  }

  int count_curr() const throw() { return count_curr_; }
  int count_max() const throw() { return count_max_; }

  void set_value (int count_max) throw ()
  {
    count_max_ = count_max;
  }

private:
  int count_max_;
  int count_curr_;
};

#endif /* CONFIG_USE_CHARM */

#endif /* CHARM_COUNTER_HPP */

