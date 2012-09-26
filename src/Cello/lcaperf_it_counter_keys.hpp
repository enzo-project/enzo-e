// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_ItCounterKeys.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-06-29
/// @brief    [\ref Lcaperf] Declaration of the ItCounterKeys class

#ifndef LCAPERF_IT_COUNTER_KEYS_HPP
#define LCAPERF_IT_COUNTER_KEYS_HPP

class ItCounterKeys {

  /// @class    ItCounterKeys
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Iterator over keys in a Counters object

public: // interface

  /// Create an ItCounterKeys object
  ItCounterKeys (Counters * counters) throw ();

  /// Delete the ItCounterKeys object
  ~ItCounterKeys () throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    p | * counters_;
    WARNING("ItCounterKeys::pup","Not pup'ing iter_");
    //    p |  iter_;
    WARNING("ItCounterKeys::pup","Not pup'ing value_");
    //    p |  value_;
  }
#endif
  
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
