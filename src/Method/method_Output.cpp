// $Id: method_Output.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Output.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Output class

#include "method.hpp"

//----------------------------------------------------------------------

Output::Output() throw ()
{
  INCOMPLETE("Output::Output");
}

//----------------------------------------------------------------------

Output::~Output() throw ()
{
  INCOMPLETE("Output::!Output");
}

//----------------------------------------------------------------------

Output::Output(const Output & output) throw ()
/// @param     output  Object being copied
{
  INCOMPLETE("Output::Output(Output)");
}

//----------------------------------------------------------------------

Output & Output::operator= (const Output & output) throw ()
/// @param     output  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("Output::operator=");
  return *this;
}

//======================================================================

void Output::write_image (Mesh * mesh, std::string filename) throw()
{
}

//----------------------------------------------------------------------

void Output::write_image (Patch * patch, bool top_level) throw()
{
}

//----------------------------------------------------------------------


void Output::write_image (Block * block, bool top_level) throw()
{
}

//----------------------------------------------------------------------

void Output::write_data (Mesh * mesh, bool top_level) throw()
{
}

//----------------------------------------------------------------------

void Output::write_data (Patch * patch, bool top_level) throw()
{
}

//----------------------------------------------------------------------

void Output::write_data (Block * block, bool top_level) throw()
{
}

//----------------------------------------------------------------------

bool Output::is_active_cycle_(int cycle) throw()
{
}

//----------------------------------------------------------------------

/// To be called at time "time" before timestep "dt" taken.  Updates
/// dt if needed.  returns (time < t_dump) && (t_dump <= time+dt) for
/// some dump time t_dump adjusts dt so that t_dump == time+dt

bool Output::is_active_time_(double time, double * p_dt) throw()
{
  bool is_active = false;
  if (time_range_.size() != 0) {
    ASSERT("Output::is_active_time_",
	   "time_range_ size must be divisible by 3",
	   time_range_.size() % 3 == 0);
    // Loop through [start,stop,step] triplets
    
    for (int i=0; i<time_range_.size(); i+=3) {
      // if triplet is a candidate, determine if time and time+dt
      // stradles a dump 
      //      @@@ LOGIC (time < t_dump) && (t_Dump <= time+dt)
      if (time < time_range_[i] && time_range_[i+2] <= time + (*p_dt)) {

	{
	  is_active = true;
	  break;
	}
	
      }
    }
    return false;
  }
}

//----------------------------------------------------------------------

