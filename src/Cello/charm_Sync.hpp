// See LICENSE_CELLO file for license and copyright information

/// @file     charm_Sync.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-12
/// @brief    [\ref Parallel] Declaration of the CHARM Sync class
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
/// First A creates and initializes a Sync variable, set to the number
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

#ifndef CHARM_SYNC_HPP
#define CHARM_SYNC_HPP

class Sync {

  /// @class    Sync
  /// @ingroup  Charm
  /// @brief    [\ref Parallel] 

 public:
   /// Create a CHARM++ "Sync" object
   Sync (int index_stop = 0) throw()
    : index_stop_(index_stop),
      index_curr_(0),
      callback_()
  {}

  /// CHARM++ pack / unpack
  void pup(PUP::er &p)
  {
    TRACEPUP;
    p | index_stop_;
    p | index_curr_;
    p | callback_;
  }

  /// Increment counter and return whether the CHARM++ parallel "sync" is done.
  inline bool next () throw()
  {
    if (index_stop_ > 0) {
      index_curr_ = (index_stop_ + (index_curr_-1) + 1) % index_stop_ + 1;  
    }
    if (index_curr_ == index_stop_) {
      index_curr_ = 0;
      return true;
    } else {
      return false;
    }
  }

  /// Return whether the Sync counter has reached the stopping value
  inline bool is_done () const throw()
  { return (index_curr_ == index_stop_);  }

  /// Set the stopping value for the counter
  inline void set_stop (int stop) throw ()
  { index_stop_ = stop; }

  /// Return the currently-set stopping value for the counter
  inline int stop () const throw ()
  { return index_stop_; }

  /// Reset the counter to 0
  inline void reset () throw () 
  { index_curr_ = 0; }

  /// Decrement the stopping value by 
  inline int operator -- () 
  { --index_stop_; return index_stop_; }

  /// Increment the stopping value by one
  inline int operator ++ () 
  { ++index_stop_; return index_stop_; }

private:

  /// Last value of the parallel sync index
  int index_stop_;

  /// Current value of the parallel sync index
  int index_curr_;

  /// Callback function
  CkCallback callback_;

};

#endif /* CHARM_SYNC_HPP */

