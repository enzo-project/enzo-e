// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleInterval.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the ScheduleInterval class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

ScheduleInterval::ScheduleInterval () throw()
  : Schedule()
{
}

//----------------------------------------------------------------------

void ScheduleInterval::set_cycle_interval
(
 int cycle_start,
 int cycle_step,
 int cycle_stop
) throw()
{
  active_ = true;
  type_ = schedule_type_cycle;

  cycle_start_ = cycle_start;
  cycle_step_  = cycle_step;
  cycle_stop_  = cycle_stop;
}

//----------------------------------------------------------------------

void ScheduleInterval::set_time_interval
  (double time_start, double time_step, double time_stop) throw()
{
  active_ = true;
  type_ = schedule_type_time;

  time_start_ = time_start;
  time_step_  = time_step;
  time_stop_  = time_stop;
}

//----------------------------------------------------------------------

bool ScheduleInterval::write_this_cycle ( int cycle, double time ) throw()
{
  bool result = false;

  if (! active_) return false;

  const double tol = 2*cello::machine_epsilon(precision_single);

  switch (type_) {

  case schedule_type_time:
    {
      bool in_range = (time_start_ <= time && time <= time_stop_);

      bool below_tol = (cello::err_abs(time_next(), time) < tol);

      result = in_range && below_tol;

    }

    break;

  case schedule_type_cycle:
    {
      bool in_range = (cycle_start_ <= cycle && cycle <= cycle_stop_);
      result = in_range && ( ((cycle-cycle_start_) % cycle_step_) == 0);

    }
    break;
  default:
    WARNING("ScheduleInterval::write_next_cycle",
	    "Unknown schedule type for active ScheduleInterval object");
  }

  return result;
}

//----------------------------------------------------------------------

double ScheduleInterval::update_timestep ( double time, double dt) const throw()
{
  if (! active_) return dt;

  double new_dt = dt;

  switch (type_) {

  case schedule_type_time:
    // time_interval_
    {

      bool in_range = (time_start_ < time && time < time_stop_);

      if (in_range) {

	double time_next = this->time_next();

	if (time < time_next && time_next < time + dt) {
	  new_dt = time_next - time;
	}
      }
    }
    break;
  default:
    break;
  }
  return new_dt;
}

