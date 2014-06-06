// See LICENSE_CELLO file for license and copyright information

/// @file     io_Schedule.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the Schedule component

#ifndef IO_SCHEDULE_HPP
#define IO_SCHEDULE_HPP

//----------------------------------------------------------------------
/// @enum     schedule_enum
/// @brief    Scheduling types

enum schedule_enum {
  schedule_type_unknown,
  schedule_type_cycle,
  schedule_type_time
};

typedef int schedule_type;

class Schedule : public PUP::able {

  /// @class    Schedule
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for various types of output

public: // functions

  /// Create an uninitialized Schedule object with the given file_name format
  Schedule() throw();

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Schedule);

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    p | active_;
    p | type_;
    p | last_;
  }

  /// Set whether the Schedule object is active or not
  void set_active(bool active) throw()
  { active_ = active; };

  /// Return whether the output object is active
  bool is_active() const throw()
  { return active_; };

  /// Set whether the Schedule object is type or not
  void set_type(int type) throw()
  {
    type_ = type; 
  };

  /// Return whether the output object is type
  int type() const throw()  { return type_; };

  /// Advance to next scheduled output
  void next () throw () { ++last_; };

public: // virtual functions

  /// Reduce timestep if next write time is between time and time + dt
  virtual double update_timestep(double time, double dt) const throw() = 0;

  /// Whether to perform IO this cycle
  virtual bool write_this_cycle ( int cycle, double time ) throw() = 0;

  virtual double time_next() const throw() = 0;

public: // static functions

  static Schedule * create (std::string var, std::string type,
			    double start,double stop,double step,
			    std::vector<double> list);

protected: // attributes

  /// Whether Schedule is currently active
  bool active_;

  schedule_type type_;

  int last_;
};

#endif /* IO_SCHEDULE_HPP */
