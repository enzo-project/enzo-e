// $Id: lcaperf_counters_deriv.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE file for license and copyright information

#ifndef LCAPERF_COUNTERS_DERIV_HPP
#define LCAPERF_COUNTERS_DERIV_HPP

/// @file     lcaperf_counters_deriv.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon May 23 11:12:54 PDT 2011
/// @brief    [\ref lcaperf] Declaration of the CountersDeriv class

class CountersDeriv : public CountersUser {

  /// @class    CountersDeriv
  /// @ingroup  lcaperf
  /// @brief    [\ref lcaperf] Derived counters

public: // interface

  /// Constructor
  CountersDeriv() throw();

  /// Destructor
  ~CountersDeriv() throw();

  /// Copy constructor
  CountersDeriv(const CountersDeriv & counters) throw();

  /// Assignment operator
  CountersDeriv & operator= (const CountersDeriv & counters) throw();

private: // attributes

};

#endif /* LCAPERF_COUNTERS_DERIV_HPP */

