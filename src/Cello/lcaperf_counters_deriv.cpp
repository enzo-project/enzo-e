// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersDeriv.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-05-27
/// @brief    Implementation of the CountersDeriv class

#include "lcaperf.hpp"

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
