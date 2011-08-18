// See LICENSE_CELLO file for license and copyright information

/// @file     field_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      time interval scheduling gets confused if multiple outputs scheduled per timestep
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Schedule class

#include "cello.hpp"

#include "field.hpp"

//----------------------------------------------------------------------

Schedule::Schedule () throw()
  : active_(false),
    schedule_type_(schedule_type_unknown),
    cycle_interval_(),
    cycle_list_(),
    time_interval_(),
    time_list_()
{
}

//----------------------------------------------------------------------

void Schedule::set_cycle_interval
(
 int cycle_start,
 int cycle_step,
 int cycle_stop
) throw()
{
  if (schedule_type_ != schedule_type_unknown) {
    WARNING ("Schedule::set_cycle_interval",
	     "Resetting Schedule scheduling");
  }

  schedule_type_ = schedule_type_cycle_interval;

  cycle_interval_.clear();

  cycle_interval_.push_back(cycle_start);
  cycle_interval_.push_back(cycle_step);
  cycle_interval_.push_back(cycle_stop);

  active_ = true;
}

//----------------------------------------------------------------------

void Schedule::set_cycle_list (std::vector<int> cycle_list) throw()
{
  if (schedule_type_ != schedule_type_unknown) {
    WARNING ("Schedule::set_cycle_list",
	     "Resetting Output scheduling");
  }

  schedule_type_ = schedule_type_cycle_list;

  cycle_list_.clear();

  cycle_list_ = cycle_list;

  active_ = true;

}

//----------------------------------------------------------------------

void Schedule::set_time_interval
  (double time_start, double time_step, double time_stop) throw()
{
  if (schedule_type_ != schedule_type_unknown) {
    WARNING ("Schedule::set_time_interval",
	     "Resetting Output scheduling");
  }

  schedule_type_ = schedule_type_time_interval;

  time_interval_.clear();

  time_interval_.push_back(time_start);
  time_interval_.push_back(time_step);
  time_interval_.push_back(time_stop);

  active_ = true;
}

//----------------------------------------------------------------------

void Schedule::set_time_list (std::vector<double> time_list) throw()
{
  if (schedule_type_ != schedule_type_unknown) {
    WARNING ("Schedule::set_time_list",
	     "Resetting Output scheduling");
  }

  schedule_type_ = schedule_type_time_list;

  time_list_.clear();

  time_list_ = time_list;

  active_ = true;
}

