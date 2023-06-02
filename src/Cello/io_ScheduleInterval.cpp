// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleInterval.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the ScheduleInterval class

#include "cello.hpp"

#include "io.hpp"

// #define DEBUG_SCHEDULE

//----------------------------------------------------------------------

ScheduleInterval::ScheduleInterval () throw()
  : Schedule(),
    cycle_start_(0),
    cycle_step_(0),
    cycle_stop_(0),
    time_start_(0.0),
    time_step_(0.0),
    time_stop_(0.0),
    seconds_start_(0.0),
    seconds_step_(0.0),
    seconds_stop_(0.0)
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

void ScheduleInterval::set_seconds_interval
  (double seconds_start, double seconds_step, double seconds_stop) throw()
{
  active_ = true;
  type_ = schedule_type_seconds;

  seconds_start_ = seconds_start;
  seconds_step_  = seconds_step;
  seconds_stop_  = seconds_stop;
}

//----------------------------------------------------------------------

namespace { // functions in anonymous namespace are local to this file

int steps_to_time_or_before_(double time, double time_start, double time_step) {
  // implicit assumption that time_start < time and time_step > 0
  int nsteps = int(std::trunc((time - time_start) / time_step));
  if ((nsteps * time_step + time_start) > time){
    return nsteps - 1;
  } else {
    return nsteps;
  }
}

};

//----------------------------------------------------------------------

bool ScheduleInterval::write_this_cycle ( int cycle, double time) throw()
{
#ifdef DEBUG_SCHEDULE
  CkPrintf ("ScheduleInterval::write_this_cycle cycle %d time %g",
	    cycle,time);
  this->print();
#endif
  double seconds = timer_.value();

  bool result = false;

  if (! active_) return false;

  const double tol = 2*cello::machine_epsilon(precision_single);

  switch (type_) {

  case schedule_type_time:
    {
      const bool in_range  = (time_start_ <= time && time <= time_stop_);

      if (in_range >= 0){
        int nsteps = steps_to_time_or_before_(time, time_start_, time_step_);
        double min_err_abs =
          std::min(cello::err_abs(nsteps * time_step_ + time_start_, time),
                   cello::err_abs((nsteps+1) * time_step_ + time_start_, time));
        result = min_err_abs < tol;
      } else {
        result = false;
      }
      // the above if-statement is a crude workaround to make this method
      // return the correct answer in the scenario when last_ isn't up-to-date.
      // It replaced the following code:
      //   const bool below_tol = (cello::err_abs(time_next(), time) < tol);
      //   result = in_range && below_tol;

    }

    break;

  case schedule_type_seconds:
    {
      const bool in_range  = (seconds_start_ <= seconds && seconds <= seconds_stop_);
      const bool past_next = (seconds > seconds_next());

      result = in_range && past_next;

    }

    break;

  case schedule_type_cycle:
    {
      const bool in_range = (cycle_start_ <= cycle && cycle <= cycle_stop_);
      const bool write_cycle = cycle_step_ ?
	( ((cycle-cycle_start_) % cycle_step_) == 0) : false;

      result = in_range && write_cycle;

    }
    break;
  default:
    WARNING("ScheduleInterval::write_next_cycle",
	    "Unknown schedule type for active ScheduleInterval object");
  }

#ifdef DEBUG_SCHEDULE
  CkPrintf ("ScheduleInterval::write_this_cycle result = %d\n",result);
#endif
  return result;
}

//----------------------------------------------------------------------

double ScheduleInterval::update_timestep ( double time, double dt)
  const throw()
{
  if (! active_) return dt;

  double new_dt = dt;

  switch (type_) {

  case schedule_type_time:
    // time_interval_
    {

      bool in_range = (time_start_ <= time && time < time_stop_);

      if (in_range) {

        int nsteps = 1 + steps_to_time_or_before_(time, time_start_,
                                                  time_step_);
        double time_next = nsteps * time_step_ + time_start_;
        // the abrove 3 lines implement a crude-workaround to the following to
        // ensure correct behavior when last_ isn't properly updated...
	//double time_next = this->time_next();

	if (time < time_next && time_next < (time + dt)) {
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

