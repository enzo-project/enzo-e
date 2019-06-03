// See LICENSE_CELLO file for license and copyright information

/// @file     cello_Sync.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-12
/// @brief    [\ref Parallel] Declaration of the Sync class
///
/// This class implements a basic count and is used to simplify
/// control flow of CHARM++ programs.  A key feature is counting
/// can begin before the stopping value is set.

#ifndef CELLO_SYNC_HPP
#define CELLO_SYNC_HPP

class Sync {

  /// @class    Sync
  /// @ingroup  Cello
  /// @brief    [\ref Parallel] 

public:
  /// Create a Sync counter object
  Sync (int index_stop = 0);

  /// CHARM++ pack / unpack
  void pup(PUP::er &p);

  /// Increment counter and return whether the counter reached the stopping value
  bool next () throw();
  
  /// Return whether the Sync counter has reached the stopping value
  bool is_done () const throw();
  
  /// Set the stopping value for the counter
  void set_stop (int stop) throw ();

  /// Increment the stopping value for the counter
  void inc_stop (int increment) throw ();

  /// Return the current counter index
  int value () const;

  /// Return the currently-set stopping value for the counter
  int stop () const throw ();

  /// Reset the counter to 0
  void reset () throw () ;
  
  /// Decrement the stopping value by one
  inline int operator -- () 
  { --index_stop_; return index_stop_; }

  /// Increment the stopping value by one
  inline int operator ++ () 
  { ++index_stop_; return index_stop_; }

  /// Decrement the stopping value by count
  inline int operator -= (int count) 
  { index_stop_ -= count; return index_stop_; }

  /// Increment the stopping value by count
  inline int operator += (int count) 
  { index_stop_ += count; return index_stop_; }

private:

  int is_done_;
  
  /// Last value of the parallel sync index
  int index_stop_;

  /// Current value of the parallel sync index
  int index_curr_;
};

#endif /* CELLO_SYNC_HPP */

