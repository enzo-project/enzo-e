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
    Schedule::pup(p);
    // NOTE: change this function whenever attributes change
    p | cycle_start_;
    p | cycle_step_;
    p | cycle_stop_;
    p | time_start_;
    p | time_step_;
    p | time_stop_;
  }

  /// Set cycle interval (start, step, stop)
  void set_cycle_interval
  (int cycle_start, int cycle_step, int cycle_stop) throw();

  /// Set time interval (start, step, stop)
  void set_time_interval
  (double time_start, double time_step, double time_stop) throw();

  /// Set seconds interval (start, step, stop)
  void set_seconds_interval
  (double seconds_start, double seconds_step, double seconds_stop) throw();

public:  // virtual functions

  /// Reduce timestep if next write time is between time and time + dt
  virtual double update_timestep(double time, double dt)  const throw();

  /// Whether to perform IO this cycle
  virtual bool write_this_cycle ( int cycle, double time) throw();

  virtual double time_next() const throw()
  { return time_start_ + (last_ + 1)*time_step_; };

  virtual double seconds_next() const throw()
  { return seconds_start_ + (last_ + 1)*seconds_step_; };

protected: // attributes

  int cycle_start_;
  int cycle_step_;
  int cycle_stop_;

  double time_start_;
  double time_step_;
  double time_stop_;

  double seconds_start_;
  double seconds_step_;
  double seconds_stop_;
};

#endif /* IO_SCHEDULE_HPP */
