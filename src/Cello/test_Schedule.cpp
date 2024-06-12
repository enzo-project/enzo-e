// See LICENSE_CELLO file for license and copyright information

/// @file     test_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Schedule class

#include "main.hpp"
#include "test.hpp"

#include "io.hpp"

#include <sys/wait.h>

static void sleep_for_2_sec() {
  // the compiler may produce warnings if we don't properly handle the exit code
  int status = system("sleep 2");

  if (status == -1) {
    ERROR("sleep_for_2_sec", "THERE WAS SOME KIND OF PROBLEM");
  } else if (WIFEXITED(status) == 0) {
    ERROR("sleep_for_2_sec", "The sleep call didn't terminate properly");
  } else if (WEXITSTATUS(status) != 0) {
    ERROR("sleep_for_2_sec",
          "The sleep command encountered some kind of error: %d.");
  }
}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ScheduleInterval");

  {
    ScheduleInterval * schedule = new ScheduleInterval;

    unit_assert (schedule != NULL);

    //--------------------------------------------------


    schedule->set_active(true);
    unit_assert(schedule->is_active());

    schedule->set_active(false);
    unit_assert(! schedule->is_active());

    schedule->set_active(true);
    unit_assert(schedule->is_active());


    unit_func("write_this_cycle");
    schedule->set_time_interval(10.0,2.0,20.0);

    unit_assert(schedule->write_this_cycle (0,  8.0) == false);
    unit_assert(schedule->write_this_cycle (0, 10.0) == true);
    unit_assert(schedule->write_this_cycle (0, 12.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 10.0) == false);
    unit_assert(schedule->write_this_cycle (0, 12.0) == true);
    unit_assert(schedule->write_this_cycle (0, 14.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 12.0) == false);
    unit_assert(schedule->write_this_cycle (0, 14.0) == true);
    unit_assert(schedule->write_this_cycle (0, 16.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 14.0) == false);
    unit_assert(schedule->write_this_cycle (0, 16.0) == true);
    unit_assert(schedule->write_this_cycle (0, 18.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 16.0) == false);
    unit_assert(schedule->write_this_cycle (0, 18.0) == true);
    unit_assert(schedule->write_this_cycle (0, 20.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 18.0) == false);
    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 22.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 22.0) == false);

    unit_func("update_timestep");

    unit_assert(schedule->update_timestep(11.0, 0.1) == 0.1);
    unit_assert(schedule->update_timestep(12.0, 0.1) == 0.1);

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

    delete schedule;
  }


  {
    ScheduleInterval * schedule = new ScheduleInterval;

    unit_assert (schedule != NULL);

    //--------------------------------------------------


    schedule->set_active(true);
    unit_assert(schedule->is_active());

    schedule->set_active(false);
    unit_assert(! schedule->is_active());

    schedule->set_active(true);
    unit_assert(schedule->is_active());


    unit_func("write_this_cycle");
    schedule->set_seconds_interval(1.0,3.0,5.0);

    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

    sleep_for_2_sec();

    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

    sleep_for_2_sec();

    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

  }
  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

