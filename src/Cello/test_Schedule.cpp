// See LICENSE_CELLO file for license and copyright information

/// @file     test_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Schedule class

#include "main.hpp"
#include "test.hpp"

#include "io.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Schedule");

  Schedule * schedule = new Schedule;

  unit_assert (schedule != NULL);

  //--------------------------------------------------


  schedule->set_active(true);
  unit_assert(schedule->is_active());
  schedule->set_active(false);
  unit_assert(! schedule->is_active());
  schedule->set_active(true);
  unit_assert(schedule->is_active());

  {
    unit_func("set_cycle_interval");

    schedule->set_cycle_interval(10,2,20);
    unit_assert(schedule->write_this_cycle (8,0.0) == false);
    unit_assert(schedule->write_this_cycle (9,0.0) == false);
    unit_assert(schedule->write_this_cycle (10,0.0) == true);
    unit_assert(schedule->write_this_cycle (11,0.0) == false);
    unit_assert(schedule->write_this_cycle (12,10.0) == true);
    unit_assert(schedule->write_this_cycle (18,0.0) == true);
    unit_assert(schedule->write_this_cycle (19,0.0) == false);
    unit_assert(schedule->write_this_cycle (20,0.0) == true);
    unit_assert(schedule->write_this_cycle (21,0.0) == false);
  }    

  unit_func("set_cycle_list");

  unit_assert (unit_incomplete);
  // std::vector<int> list;
  // list.push_back(12);
  // list.push_back(-18);
  // list.push_back(7);
  // schedule->set_cycle_list(list);

  unit_func("set_time_interval");

  unit_assert (unit_incomplete);

  // schedule->set_time_interval(13.3, 2.0, 19);


  unit_func("set_time_list");

  unit_assert (unit_incomplete);
  // std::vector<double> list_double;
  // list.push_back(7/3);
  // list.push_back(-18.9);
  // list.push_back(7);
  // schedule->set_time_list (list_double);


  unit_func("this_cycle ");

  unit_assert (unit_incomplete);

  unit_func("t = schedule->update_timestep");

  unit_assert (unit_incomplete);

  unit_func("set_skip_cycle ");
  unit_assert (unit_incomplete);
  unit_func("set_skip_time ");
  unit_assert (unit_incomplete);


  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

