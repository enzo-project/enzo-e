// $Id: lcaperf_it_counter_keys.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE file for license and copyright information

#ifndef LCAPERF_IT_COUNTER_KEYS_HPP
#define LCAPERF_IT_COUNTER_KEYS_HPP

/// @file     lcaperf_it_counter_keys.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jun 29 18:31:41 PDT 2011
/// @brief    [\ref lcaperf] Declaration of the ItCounterKeys class

class ItCounterKeys {

  /// @class    ItCounterKeys
  /// @ingroup  lcaperf
  /// @brief    [\ref lcaperf] Iterator over keys in a Counters object

public: // interface

  /// Create an ItCounterKeys object
  ItCounterKeys (Counters * counters) throw ();

  /// Delete the ItCounterKeys object
  ~ItCounterKeys () throw ();
  
  /// Iterate through all local CounterKeys in the Mesh
  const char * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

  /// Return the counter value of the corresponding key
  long long * value() const throw();

private: // attributes

  /// The counters being iterated over
  Counters * counters_;
  
  /// Iterator over Counters keys
  std::map<std::string,long long *>::iterator iter_;

  /// Copy of counter value array for value()
  long long * value_;
};

#endif
