// See LICENSE_CELLO file for license and copyright information

/// @file     io_ScheduleList.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the ScheduleList component

#ifndef IO_SCHEDULELIST_HPP
#define IO_SCHEDULELIST_HPP
class ScheduleList : public Schedule {

  /// @class    ScheduleList
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for various types of output

public: // functions

  /// Create an uninitialized ScheduleList object with the given file_name format
  ScheduleList() throw();


  /// CHARM++ PUP::able declaration
  PUPable_decl(ScheduleList);

  /// CHARM++ migration constructor
  ScheduleList(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    Schedule::pup(p);
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    p | cycle_list_;
    p | time_list_;
  }

  /// Set cycle list
  void set_cycle_list (std::vector<int> cycle_list) throw();

  /// Set time list
  void set_time_list (std::vector<double> time_list) throw();

  /// Set seconds list
  void set_seconds_list (std::vector<double> seconds_list) throw();

public: // virtual functions

  /// Reduce timestep if next write time is between time and time + dt
  virtual double update_timestep(double time, double dt) const throw();

  /// Whether to perform IO this cycle
  virtual bool write_this_cycle ( int cycle, double time ) throw();

  virtual double time_next() const throw()
  { return ((int(time_list_.size()) > last_+1) ? time_list_.at(last_+1) : -1.0); }

  virtual double seconds_next() const throw()
  { return ((int(seconds_list_.size()) > last_+1) ? seconds_list_.at(last_+1) : -1.0); }

protected: // attributes
  /// List of cycles to perform schedule
  std::vector<int> cycle_list_;

  /// List of times to perform schedule
  std::vector<double> time_list_;

  /// List of seconds clock times to perform schedule
  std::vector<double> seconds_list_;

  /// Previous seconds time
  double seconds_prev_;
};

#endif /* IO_SCHEDULE_HPP */
