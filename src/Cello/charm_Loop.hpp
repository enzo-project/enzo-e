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
///
/// A::foo()
/// {
/// }
///
/// B::p_foo()
/// {
/// }
/// 
/// A::s_foo()
/// {
/// }

#ifndef CHARM_LOOP_HPP
#define CHARM_LOOP_HPP

class Loop {

  /// @class    Loop
  /// @ingroup  Charm
  /// @brief    [\ref Charm] 

 public:
   /// Create a CHARM++ "Loop" object
   Loop (int index_stop = 0
	) throw()
    : index_stop_(index_stop),
      index_curr_(0)
  {}

  /// CHARM++ pack / unpack
#ifdef CONFIG_USE_CHARM
  void pup(PUP::er &p)
  {
    TRACEPUP;
    p | index_stop_;
    p | index_curr_;
  }
#endif

  /// Return whether the CHARM++ parallel "loop" is done
  inline bool done (int index = 1) throw()
  {
    if (index_stop_ > 0) {
      index_curr_ = (index_stop_ + index_curr_ - index) % index_stop_;  
    }
    return index_curr_ == 0;
  }

  /// Return the current CHARM++ parallel "loop" index
  inline int index() const throw() { return index_curr_; }
  /// Return the upper-limit on the CHARM++ parallel "loop"
  inline int stop() const throw()  { return index_stop_; }
  /// Access to the upper-limit on the CHARM++ parallel "loop"
  inline int & stop () throw ()    { return index_stop_; }

private:

  /// Last value of the parallel loop index
  int index_stop_;

  /// Current value of the parallel loop index
  int index_curr_;

};

#endif /* CHARM_LOOP_HPP */

