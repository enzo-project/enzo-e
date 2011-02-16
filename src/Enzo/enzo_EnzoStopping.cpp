// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @todo     Remove repeated creation / deletion of CellWidth[]
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoStopping user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoStopping::EnzoStopping 
(
 Parameters * parameters,
 EnzoDescr  * enzo
)
  : Stopping(),
    enzo_(enzo),
    cycle_stop_(-1),
    time_stop_(-1.0)
{
  //--------------------------------------------------
  parameters->set_current_group ("Stopping");
  //--------------------------------------------------

  cycle_stop_ = parameters->value_integer("cycle",1000);
  time_stop_  = parameters->value_scalar("time",2.5);
}

//----------------------------------------------------------------------

void EnzoStopping::update_block (DataBlock * block) throw()
{
  INCOMPLETE("EnzoStopping::update_block","not implemented");
}

//----------------------------------------------------------------------

bool EnzoStopping::complete () throw()
{
  ASSERT("EnzoStopping::complete",
	 "Neither Stopping::time_stop nor Stopping::cycle_stop initialized",
	 time_stop_ != -1.0 || cycle_stop_ != -1);
  return 
    (enzo_->Time >= time_stop_ ) ||
    (enzo_->CycleNumber >= cycle_stop_);
}

