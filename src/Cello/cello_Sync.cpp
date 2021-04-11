// See LICENSE_CELLO file for license and copyright information

/// @file     cello_Sync.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-05-29
/// @brief    [\ref Parallel] Implementation of the Sync class
///

#include "cello.hpp"
#include "charm.hpp"

// #define TRACE_SYNC
// #define CHECK_OVER_COUNT

#ifdef TRACE_SYNC
#  undef TRACE_SYNC
#  define TRACE_SYNC(MSG)                                               \
  CkPrintf ("TRACE_SYNC %s :%d %d/%d %d\n",                             \
            MSG,__LINE__,index_curr_,index_stop_,is_done_);\
  fflush(stdout);
#else
#  define TRACE_SYNC(MSG) /* ... */
#endif

Sync::Sync (int index_stop)
  : is_done_(0),
    index_stop_(index_stop),
    index_curr_(0),
    state_(RefreshState::INACTIVE)
{}

//----------------------------------------------------------------------

void Sync::pup(PUP::er &p)
{
  TRACEPUP;
  p | is_done_;
  p | index_stop_;
  p | index_curr_;
  p | state_;
}

//----------------------------------------------------------------------

bool Sync::next () throw()
{
  advance_();
  check_done_();
  TRACE_SYNC("next");
#ifdef CHECK_OVER_COUNT
  return is_done_;
#else  
  return (index_curr_ == 0) && is_done_;
#endif  
}

//----------------------------------------------------------------------

void Sync::advance () throw()
{
  advance_();
  check_done_();
  TRACE_SYNC("advance");
}

//----------------------------------------------------------------------

void Sync::advance_() throw()
{
  ++ index_curr_;
}

//----------------------------------------------------------------------

void Sync::check_done_() throw()
{
  if (index_stop_ > 0) {
    if (index_curr_ == index_stop_) {
      // reached stopping value
#ifndef CHECK_OVER_COUNT      
      index_curr_ = 0;
#endif      
      is_done_ = true;
    }
    if (index_curr_ > index_stop_) {
      // exceded stopping value: error!
      ERROR2("Sync::check_done_()",
             "Incrementing sync counter %d beyond limit %d",
             index_curr_,index_stop_);
    }
  }
}
//----------------------------------------------------------------------

bool Sync::is_done () const throw()
{ return is_done_; }

//----------------------------------------------------------------------

void Sync::set_stop (int stop) throw ()
{
  index_stop_ = stop;
  TRACE_SYNC("set_stop");
}


//----------------------------------------------------------------------

void Sync::inc_stop (int increment) throw ()
{
  index_stop_ += increment;
  TRACE_SYNC("inc_stop");
}

//----------------------------------------------------------------------
int Sync::value () const
{ return index_curr_; }

//----------------------------------------------------------------------
int Sync::stop () const throw ()
{ return index_stop_; }

//----------------------------------------------------------------------
void Sync::reset () throw () 
{ index_curr_ = 0;
  index_stop_ = 0;
  is_done_ = 0;
  TRACE_SYNC("reset");
}

