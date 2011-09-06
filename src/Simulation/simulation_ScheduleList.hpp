// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_Schedule.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Simulation] Declaration for the ScheduleList component

#ifndef SIMULATION_SCHEDULE_HPP
#define SIMULATION_SCHEDULE_HPP

class ScheduleList : public Schedule {

  /// @class    ScheduleList
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] define interface for various types of simulation output

public: // functions

  /// Create an uninitialized ScheduleList object with the given file_name format
  ScheduleList() throw();

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

  /// Set whether the ScheduleList object is active or not
  void set_active(bool active) throw()
  { active_ = active; };

  /// Return whether the output object is active
  bool is_active() const throw()
  { return active_; };

  /// Whether to perform IO this cycle
  bool write_this_cycle ( int cycle, double time ) throw();

  /// Reduce timestep if next write time is greater than time + dt
  double update_timestep(double time, double dt) const throw();

protected: // attributes

  /// Whether ScheduleList is currently active
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

  /// Index of time or cycle interval or list for next output
  size_t index_;

};

#endif /* SIMULATION_SCHEDULE_HPP */
