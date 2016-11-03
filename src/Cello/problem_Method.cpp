// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-09-04
/// @brief    

#include "problem.hpp"

double Method::courant_global = 1.0;

//----------------------------------------------------------------------

Method::~Method() throw()
{
  delete schedule_;
  for (size_t i=0; i<refresh_list_.size(); i++) {
    delete refresh_list_[i];
    refresh_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

void Method::pup (PUP::er &p)
{ TRACEPUP;
  PUP::able::pup(p);
  p | refresh_list_;
  p | schedule_; // pupable
  p | courant_;
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule) throw()
{ 
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//======================================================================

