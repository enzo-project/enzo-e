// $Id: counters_deriv.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE file for license and copyright information

/// @file     counters_deriv.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-27
/// @brief    Implementation of the CountersDeriv class

#include "performance.hpp"

namespace lca {

//----------------------------------------------------------------------

CountersDeriv::CountersDeriv() throw ()
{
}

//----------------------------------------------------------------------

CountersDeriv::~CountersDeriv() throw ()
{
}

//----------------------------------------------------------------------

CountersDeriv::CountersDeriv(const CountersDeriv & counters) throw ()
  : CountersUser(counters)
/// @param     counters  Object being copied
{
}

//----------------------------------------------------------------------

CountersDeriv & CountersDeriv::operator= (const CountersDeriv & counters) throw ()
/// @param     counters  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================
}
