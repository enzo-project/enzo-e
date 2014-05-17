// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleInterval.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the ScheduleInterval class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

ScheduleInterval::ScheduleInterval () throw()
  : Schedule(),
    cycle_interval_(),
    time_interval_()
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
  if (type_ != schedule_type_unknown) {
    WARNING ("ScheduleInterval::set_cycle_interval",
	     "Resetting Schedule scheduling");
  }

  cycle_interval_.clear();

  cycle_interval_.push_back(cycle_start);
  cycle_interval_.push_back(cycle_step);
  cycle_interval_.push_back(cycle_stop);

  active_ = true;
}

//----------------------------------------------------------------------

void ScheduleInterval::set_time_interval
  (double time_start, double time_step, double time_stop) throw()
{
  if (type_ != schedule_type_unknown) {
    WARNING ("Schedule::set_time_interval",
	     "Resetting Output scheduling");
  }

  time_interval_.clear();

  time_interval_.push_back(time_start);
  time_interval_.push_back(time_step);
  time_interval_.push_back(time_stop);

  active_ = true;
}

//----------------------------------------------------------------------

bool ScheduleInterval::write_this_cycle ( int cycle, double time ) throw()
{
  bool result = false;

  if (! active_) return false;

  const double tol = 2*cello::machine_epsilon(precision_single);

  double time_start;
  double time_step;
  double time_stop;

  int cycle_start;
  int cycle_step;
  int cycle_stop;

  switch (type_) {

  case schedule_type_time:
    {
      time_start = time_interval_[0];
      time_step  = time_interval_[1];
      time_stop  = time_interval_[2];

      bool in_range = (time_start <= time && time <= time_stop);

      double ratio = (time - time_start) / time_step;
      bool below_tol = (cello::err_abs(round(ratio),ratio) < tol);

      result = in_range && below_tol;

    }

    break;

  case schedule_type_cycle:
    {
      cycle_start = cycle_interval_[0];
      cycle_step  = cycle_interval_[1];
      cycle_stop  = cycle_interval_[2];

      bool in_range = (cycle_start <= cycle && cycle <= cycle_stop);
      result = in_range && ( ((cycle-cycle_start) % cycle_step) == 0);

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

  const double time_next  = time + dt;

  double time_start;
  double time_step;
  double time_stop;
  double time_dump;

  switch (type_) {

  case schedule_type_time:
    // time_interval_
    {
      time_start = time_interval_[0];
      time_step  = time_interval_[1];
      time_stop  = time_interval_[2];

      double ratio      = (time - time_start) / time_step;
      double ratio_next = (time_next - time_start) / time_step;

      if ((round(ratio) == round(ratio_next)) &&
	  ratio < round(ratio) && ratio_next > round(ratio)) {
	time_dump = time_start + round(ratio)*time_step;
	new_dt = time_dump - time;
      }
    }
    break;
  default:
    break;
  }
  return new_dt;
}

