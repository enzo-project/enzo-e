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
    type_(schedule_type_unknown)
{
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
  bool var_cycle,var_time;

  var_cycle = (var == "cycle");
  var_time  = (var == "time");

  bool type_interval,type_list;

  type_interval = type == "interval";
  type_list     = type == "list";

    
  Schedule * schedule = 0;

  if (type_interval) {

    schedule = new ScheduleInterval;

    if (var_cycle) {

      ((ScheduleInterval * )schedule)->set_cycle_interval(start,step,stop);
      
    } else if (var_time) {

      ((ScheduleInterval * )schedule)->set_time_interval(start,step,stop);
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

    }

  }

  if (var_cycle) schedule->set_type(schedule_type_cycle);
  if (var_time)  schedule->set_type(schedule_type_time);

  return schedule;

}
