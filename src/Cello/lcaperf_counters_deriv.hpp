// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersDeriv.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-05-23
/// @brief    [\ref Lcaperf] Declaration of the CountersDeriv class
///

#ifndef LCAPERF_COUNTERS_DERIV_HPP
#define LCAPERF_COUNTERS_DERIV_HPP

class CountersDeriv : public CountersUser {

  /// @class    CountersDeriv
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Derived counters

public: // interface

  /// Constructor
  CountersDeriv() throw();

  /// Destructor
  ~CountersDeriv() throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    CountersUser::pup(p);
  }
#endif

  /// Copy constructor
  CountersDeriv(const CountersDeriv & counters) throw();

  /// Assignment operator
  CountersDeriv & operator= (const CountersDeriv & counters) throw();

private: // attributes

};

#endif /* LCAPERF_COUNTERS_DERIV_HPP */

