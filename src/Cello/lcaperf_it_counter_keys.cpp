// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_ItCounterKeys.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-19
/// @brief    Implementation of the ItCounterKeys abstract base class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

namespace lca {

//----------------------------------------------------------------------

ItCounterKeys::ItCounterKeys ( Counters * counters ) throw ()
  : counters_(counters),
    iter_(counters->global_.begin()),
    value_(0)
{}

//----------------------------------------------------------------------

ItCounterKeys::~ItCounterKeys ( ) throw ()
{}

//----------------------------------------------------------------------

const char * ItCounterKeys::operator++ () throw()
{
  const char * key;
  if (iter_ == counters_->global_.end()) {
    key = 0;
    value_ = 0;
    iter_ = counters_->global_.begin();
  } else {
    key    = iter_->first.c_str();
    value_ = iter_->second;
    ++iter_;
  }
  return key;
}

//----------------------------------------------------------------------

bool ItCounterKeys::done () const throw()
{
  return iter_ == counters_->global_.begin();
}

//----------------------------------------------------------------------

long long * ItCounterKeys::value () const throw()
{
  return value_;
}
}
