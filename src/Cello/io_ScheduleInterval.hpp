// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleInterval.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the ScheduleInterval component

#ifndef IO_SCHEDULEINTERVAL_HPP
#define IO_SCHEDULEINTERVAL_HPP

class ScheduleInterval : public Schedule {

  /// @class    ScheduleInterval
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for various types of output

public: // functions

  /// Create an uninitialized ScheduleInterval object with the given file_name format
  ScheduleInterval() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(ScheduleInterval);

  /// CHARM++ migration constructor
  ScheduleInterval(CkMigrateMessage *m) {}


  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    ScheduleInterval::pup(p);
    // NOTE: change this function whenever attributes change
    p | cycle_interval_;
    p | time_interval_;
  }

  /// Set cycle interval (start, step, stop)
  void set_cycle_interval
  (int cycle_start, int cycle_step, int cycle_stop) throw();

  /// Set time interval (start, step, stop)
  void set_time_interval
  (double time_start, double time_step, double time_stop) throw();

public:  // virtual functions

  /// Reduce timestep if next write time is between time and time + dt
  virtual double update_timestep(double time, double dt) const throw();

  /// Whether to perform IO this cycle
  bool write_this_cycle ( int cycle, double time ) throw();

protected: // attributes

  /// cycle start, step, and stop of schedule
  std::vector<int> cycle_interval_;

  /// time start, step, and stop of schedule
  std::vector<double> time_interval_;
};

#endif /* IO_SCHEDULE_HPP */
