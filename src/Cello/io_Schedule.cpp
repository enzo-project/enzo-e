// See LICENSE_CELLO file for license and copyright information

/// @file     io_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Schedule base class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

Schedule::Schedule () throw()
  : active_(false),
    type_(schedule_type_unknown),
    last_(-1),
    timer_()
{
  timer_.start();
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

  const bool var_cycle   = (var == "cycle");
  const bool var_time    = (var == "time");
  const bool var_seconds = (var == "seconds");
  const bool var_minutes = (var == "minutes");
  const bool var_hours   = (var == "hours");

  const bool type_interval = (type == "interval");
  const bool type_list     = (type == "list");

  Schedule * schedule = 0;

  if (type_interval) {

    schedule = new ScheduleInterval;

    if (var_cycle) {

      ((ScheduleInterval * )schedule)->set_cycle_interval(start,step,stop);

    } else if (var_time) {

      ((ScheduleInterval * )schedule)->set_time_interval(start,step,stop);

    } else if (var_seconds) {

      ((ScheduleInterval * )schedule)->
        set_seconds_interval(start,step,stop);

    } else if (var_minutes) {

      ((ScheduleInterval * )schedule)->
        set_seconds_interval(60.0*start,60.0*step,60.0*stop);

    } else if (var_hours) {

      ((ScheduleInterval * )schedule)->
        set_seconds_interval(3600.0*start,3600.0*step,3600.0*stop);

    }

  } else if (type_list) {

    schedule = new ScheduleList;

    if (var_cycle) {

      int size = list.size();
      std::vector<int> list_int;
      list_int.resize(size);
      for (int i=0; i<size; i++) {
	list_int[i] = (int) list[i];
      }

      ((ScheduleList * )schedule)->set_cycle_list(list_int);

    } else if (var_time) {

      ((ScheduleList * )schedule)->set_time_list(list);

    } else if (var_seconds) {

      ((ScheduleList * )schedule)->set_seconds_list(list);

    } else if (var_minutes) {

      for (std::size_t i=0; i<list.size(); i++) list[i] *= 60.0;

      ((ScheduleList * )schedule)->set_seconds_list(list);

    } else if (var_hours) {

      for (std::size_t i=0; i<list.size(); i++) list[i] *= 3600.0;

      ((ScheduleList * )schedule)->set_seconds_list(list);

    }

  }

  ASSERT ("Schedule::create",
	  "Bad shedule parameters",
	  schedule != NULL);

  if (var_cycle) schedule->set_type(schedule_type_cycle);
  if (var_time)  schedule->set_type(schedule_type_time);
  if (var_seconds || var_minutes || var_hours)
    schedule->set_type(schedule_type_seconds);

  return schedule;

}
