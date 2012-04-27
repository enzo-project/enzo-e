// $Id: it_counter_keys.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE file for license and copyright information

/// @file     it_counter_keys.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu May 19 2011
/// @brief    Implementation of the ItCounterKeys abstract base class

#include "performance.hpp"

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
  DEBUG1("counters_ = %p",counters_);
  const char * key;
  if (iter_ == counters_->global_.end()) {
  DEBUG0;
    key = 0;
    value_ = 0;
    iter_ = counters_->global_.begin();
  } else {
  DEBUG0;
    key    = iter_->first.c_str();
    value_ = iter_->second;
    ++iter_;
  }
  DEBUG0;
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
