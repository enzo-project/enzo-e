// See LICENSE_CELLO file for license and copyright information

/// @file     charm_Loop.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-12
/// @brief    [\ref Charm] Declaration of the CHARM Loop class
///
/// This class is used to simplify control flow of CHARM++ programs.
///
///    A::foo()       o        o        o      |
///                  /|\      /|\      /|\     |
///    B::p_foo()   o o o    o o o    o o o    |
///                  \|/      \|/      \|/     |
///    A::s_foo()     o        o        o      |
///
///
/// First A creates and initializes a Loop variable, set to the number
/// of expected "returns" into A::s_foo() from B::p_foo() (which may
/// or may not be the number of calls from A to B).

#ifndef CHARM_LOOP_HPP
#define CHARM_LOOP_HPP

#ifdef CONFIG_USE_CHARM
#include "charm++.h"
class Loop {

  /// @class    Loop
  /// @ingroup  Charm
  /// @brief    [\ref Charm] 

 public:
   /// Create a CHARM++ "Loop" object
   Loop (int index_max = 0
	) throw()
    : index_max_(index_max),
      index_curr_(0)
  {}

  void pup(PUP::er &p)
  {
    p | index_max_;
    p | index_curr_;
  }


  inline bool done (int index = 1) throw()
  {
    if (index_max_ > 0) {
      index_curr_ = (index_max_ + index_curr_ - index) % index_max_;  
    }
    return remaining() == 0;
  }

  inline int remaining() throw ()
  {
    return index_curr_;
  }

  int index_curr() const throw() { return index_curr_; }
  int index_max() const throw() { return index_max_; }

  void set_max (int index_max) throw ()
  {
    index_max_ = index_max;
  }
  void inc_max () throw ()
  {
    ++index_max_;
  }

private:
  int index_max_;
  int index_curr_;
};

#endif /* CONFIG_USE_CHARM */

#endif /* CHARM_LOOP_HPP */

