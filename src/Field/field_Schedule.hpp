// See LICENSE_CELLO file for license and copyright information

/// @file     field_Schedule.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @todo     Extract Schedule object for scheduling output
/// @brief    [\ref Field] Declaration for the Schedule component

#ifndef FIELD_SCHEDULE_HPP
#define FIELD_SCHEDULE_HPP

//----------------------------------------------------------------------
/// @enum     schedule_type
/// @brief    Scheduling types

enum schedule_type {
  schedule_type_unknown,
  schedule_type_cycle_interval,
  schedule_type_cycle_list,
  schedule_type_time_interval,
  schedule_type_time_list
};

class Schedule {

  /// @class    Schedule
  /// @ingroup  Field
  /// @brief    [\ref Field] define interface for various types of simulation output

public: // functions

  /// Create an uninitialized Schedule object with the given file_name format
  Schedule() throw();

  /// Set cycle interval (start, step, stop)
  void set_cycle_interval
  (int cycle_start, int cycle_step, int cycle_stop) throw();

  /// Set cycle list
  void set_cycle_list (std::vector<int> cycle_list) throw();

  /// Set time interval (start, step, stop)
  void set_time_interval
  (double time_start, double time_step, double time_stop) throw();

  /// Set time list
  void set_time_list (std::vector<double> time_list) throw();

  /// Set whether the Schedule object is active or not
  void set_active(bool active) throw()
  { active_ = active; };
    
  /// Return whether the output object is active
  bool is_active() const throw()
  { return active_; };

protected: // attributes

  /// Whether Schedule is currently active
  bool active_;

  /// Schedule type of the Schedule object
  schedule_type schedule_type_;

  /// cycle start, step, and stop of schedule
  std::vector<int> cycle_interval_;

  /// List of cycles to perform schedule
  std::vector<int> cycle_list_;

  /// time start, step, and stop of schedule
  std::vector<double> time_interval_;

  /// List of times to perform schedule
  std::vector<double> time_list_;
};

#endif /* FIELD_SCHEDULE_HPP */
