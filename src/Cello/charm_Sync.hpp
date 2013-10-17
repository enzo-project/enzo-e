// See LICENSE_CELLO file for license and copyright information

/// @file     charm_Sync.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-12
/// @brief    [\ref Charm] Declaration of the CHARM Sync class
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
  /// @brief    [\ref Charm] 

 public:
   /// Create a CHARM++ "Sync" object
   Sync (int index_stop = 0
	) throw()
    : index_stop_(index_stop),
      index_curr_(0)
  {}

  /// CHARM++ pack / unpack
  void pup(PUP::er &p)
  {
    TRACEPUP;
    p | index_stop_;
    p | index_curr_;
  }

  /// Increment counter and return whether the CHARM++ parallel "sync" is done.
  inline bool next (int index = 1) throw()
  {
    if (index_stop_ > 0) {
      index_curr_ = (index_stop_ + (index_curr_-1) + index) % index_stop_ + 1;  
    }
    return index_curr_ == index_stop_;
  }

  inline bool is_done () const throw()
  {
    return (index_curr_ == index_stop_);
  }

  inline int operator = (int value) { index_stop_ = value; return index_stop_; }
  inline int operator -= (int value) { index_stop_ -= value; return index_stop_; }
  inline int operator += (int value) { index_stop_ += value; return index_stop_; }
  inline int operator -- () { --index_stop_; return index_stop_; }
  inline int operator ++ () { ++index_stop_; return index_stop_; }

  /// Return the current CHARM++ parallel "sync" index
  inline int index() const throw() { return index_curr_; }
  /// Set the current index
  void set_index(int value) { index_curr_ = value; }
  void add_index(int value) { index_curr_ += value; }
  /// Return the upper-limit on the CHARM++ parallel "sync"
  inline int stop() const throw()  { return index_stop_; }
  /// Access to the upper-limit on the CHARM++ parallel "sync"
  inline void set_stop (int stop) throw ()    { index_stop_ = stop; }
  inline void add_stop (int inc = 1) throw () { index_stop_ += inc; }
  inline void clear () throw () { index_curr_ = 0; }

private:

  /// Last value of the parallel sync index
  int index_stop_;

  /// Current value of the parallel sync index
  int index_curr_;

};

#endif /* CHARM_SYNC_HPP */

