// See LICENSE_CELLO file for license and copyright information

/// @file     cello_Sync.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-05-29
/// @brief    [\ref Parallel] Implementation of the Sync class
///

#include "cello.hpp"
#include "charm.hpp"

Sync::Sync (int index_stop)
  : is_done_(0),
    index_stop_(index_stop),
    index_curr_(0)
{}

//----------------------------------------------------------------------

void Sync::pup(PUP::er &p)
{
  TRACEPUP;
  p | is_done_;
  p | index_stop_;
  p | index_curr_;
}

//----------------------------------------------------------------------

bool Sync::next () throw()
{
  advance_();
  check_done_();
  return (index_curr_ == 0 && is_done_ == true);
}


//----------------------------------------------------------------------

void Sync::advance () throw()
{
  advance_();
  check_done_();
}

//----------------------------------------------------------------------

void Sync::advance_() throw()
{
  if (index_stop_ > 0) {
    index_curr_ = (index_stop_ + (index_curr_-1) + 1) % index_stop_ + 1;  
  } else {
    // stop is not known yet
    ++ index_curr_;
  }
}

//----------------------------------------------------------------------

void Sync::check_done_() throw()
{
  if ( (index_curr_ == index_stop_) && 
       (index_stop_ > 0) ) {
    index_curr_ = 0;
    is_done_ = true;
  }
}
//----------------------------------------------------------------------

bool Sync::is_done () const throw()
{ return is_done_;  }

//----------------------------------------------------------------------

void Sync::set_stop (int stop) throw ()
{ index_stop_ = stop; }


//----------------------------------------------------------------------

void Sync::inc_stop (int increment) throw ()
{ index_stop_ += increment; }

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
  is_done_ = 0;}

