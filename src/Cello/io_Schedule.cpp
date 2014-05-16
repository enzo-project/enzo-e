// See LICENSE_CELLO file for license and copyright information

/// @file     io_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Schedule class

#include "cello.hpp"

#include "io.hpp"

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

//----------------------------------------------------------------------

bool Schedule::write_this_cycle ( int cycle, double time ) throw()
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

  // Check if this cycle or time is to be skipped

  bool skip = false;

  switch (schedule_type_) {
  case schedule_type_time_list:
  case schedule_type_time_interval:
    for (size_t i=0; i<time_skip_.size(); i++) {
      if (time_skip_.at(i)==time) skip = true;
    }

    break;
  case schedule_type_cycle_list:
  case schedule_type_cycle_interval:
    for (size_t i=0; i<cycle_skip_.size(); i++) {
      if (cycle_skip_.at(i)==cycle) skip = true;
    }
    break;
  default:
    WARNING("Schedule::write_next_cycle",
	    "Unknown schedule type for active Schedule object");
  }

  if (skip) {
    return false;
  }

  switch (schedule_type_) {

  case schedule_type_time_list:

    // time_list_
    for (size_t i=0; i<time_list_.size(); i++) {
      if (cello::err_abs(time,time_list_[i]) < tol) {
	result = true;
	break;
      }
    }
    break;

  case schedule_type_time_interval:
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

  case schedule_type_cycle_list:
    for (size_t i=0; i<cycle_list_.size(); i++) {
      if (cycle == cycle_list_[i]) {
	result = true;
	break;
      }
    }
    break;

  case schedule_type_cycle_interval:
    {
      cycle_start = cycle_interval_[0];
      cycle_step  = cycle_interval_[1];
      cycle_stop  = cycle_interval_[2];

      bool in_range = (cycle_start <= cycle && cycle <= cycle_stop);
      result = in_range && ( ((cycle-cycle_start) % cycle_step) == 0);

    }
    break;
  default:
    WARNING("Schedule::write_next_cycle",
	    "Unknown schedule type for active Schedule object");
  }

  return result;
}

//----------------------------------------------------------------------

double Schedule::update_timestep ( double time, double dt) const throw()
{
  if (! active_) return dt;

  double new_dt = dt;

  const double time_next  = time + dt;

  double time_start;
  double time_step;
  double time_stop;
  double time_dump;

  switch (schedule_type_) {

  case schedule_type_time_list:

    for (size_t i=0; i<cycle_list_.size(); i++) {
      time_dump = time_list_[i];
      if (time < time_dump && time_dump < time_next) {
	new_dt = time_dump - time;
	break;
      }
    }
    break;

  case schedule_type_time_interval:
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

//======================================================================

// static
Schedule * Schedule::create 
(std::string var,
 std::string type,
 double start,
 double stop,
 double step,
 std::vector<double> list)
{
  Schedule * schedule = new Schedule;

  bool var_cycle,var_time;

  var_cycle = (var == "cycle");
  var_time  = (var == "time");

  bool type_interval,type_list;

  type_interval = type == "interval";
  type_list     = type == "list";

    
  if (type_interval) {

    if (var_cycle) {

      schedule->set_cycle_interval(start,step,stop);
      
    } else if (var_time) {

      schedule->set_time_interval(start,step,stop);
    }

  } else if (type_list) {

    if (var_cycle) {

      int size = list.size();
      std::vector<int> list_int;
      list_int.resize(size);
      for (int i=0; i<size; i++) {
	list_int[i] = (int) list[i];
      }
      schedule->set_cycle_list(list_int);

    } else if (var_time) {

      schedule->set_time_list(list);

    }

  }
  
  return schedule;

}
