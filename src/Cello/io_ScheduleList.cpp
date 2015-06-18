// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleList.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the ScheduleList class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

ScheduleList::ScheduleList () throw()
  : Schedule(),
    cycle_list_(),
    time_list_()
{
}

//----------------------------------------------------------------------

void ScheduleList::set_cycle_list (std::vector<int> cycle_list) throw()
{
  if (type_ != schedule_type_unknown) {
    WARNING ("ScheduleList::set_cycle_list",
	     "Resetting Output scheduling");
  }

  type_ = schedule_type_cycle;

  cycle_list_.clear();

  cycle_list_ = cycle_list;

  active_ = true;

}

//----------------------------------------------------------------------

void ScheduleList::set_time_list (std::vector<double> time_list) throw()
{

  if (type_ != schedule_type_unknown) {
    WARNING ("ScheduleList::set_time_list",
	     "Resetting Output scheduling");
  }

  type_ = schedule_type_time;

  time_list_.clear();

  time_list_ = time_list;

  active_ = true;
}

//----------------------------------------------------------------------

void ScheduleList::set_seconds_list (std::vector<double> seconds_list) throw()
{

  if (type_ != schedule_type_unknown) {
    WARNING ("ScheduleList::set_seconds_list",
	     "Resetting Output scheduling");
  }

  type_ = schedule_type_seconds;

  seconds_list_.clear();

  seconds_list_ = seconds_list;

  active_ = true;
}

//----------------------------------------------------------------------

bool ScheduleList::write_this_cycle ( int cycle, double time ) throw()
{
  bool result = false;

  if (! active_) return false;

  const double tol = 2*cello::machine_epsilon(precision_single);

  switch (type_) {

  case schedule_type_time:


    // time_list_
    for (size_t i=0; i<time_list_.size(); i++) {
      if (cello::err_abs(time,time_list_[i]) < tol) {
	result = true;
	break;
      }
    }
    break;

  case schedule_type_cycle:
    for (size_t i=0; i<cycle_list_.size(); i++) {
      if (cycle == cycle_list_[i]) {
	result = true;
	break;
      }
    }
    break;

  default:
    WARNING("ScheduleList::write_next_cycle",
	    "Unknown schedule type for active ScheduleList object");
  }

  return result;
}

//----------------------------------------------------------------------

double ScheduleList::update_timestep ( double time, double dt) const throw()
{
  if (! active_) return dt;

  double new_dt = dt;

  const double time_next  = time + dt;

  double time_dump;

  switch (type_) {

  case schedule_type_time:

    for (size_t i=0; i<cycle_list_.size(); i++) {
      time_dump = time_list_[i];
      if (time < time_dump && time_dump < time_next) {
	new_dt = time_dump - time;
	break;
      }
    }
    break;
  default:
    break;
  }
  return new_dt;
}

