// See LICENSE_CELLO file for license and copyright information

/// @file     data_State.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2023-12-29
/// @brief    [\ref Data] Declaration of the State class

#ifndef DATA_STATE_HPP
#define DATA_STATE_HPP

class State {

  /// @class    State
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  State() throw()
  : cycle_(0),
    time_(0.0),
    dt_(0.0),
    stopping_(false),
    method_time_(),
    method_dt_()
  {
  }

  /// Constructor
  State(int cycle, double time, double dt, bool stopping) throw()
    : cycle_(cycle),
      time_(time),
      dt_(dt),
      stopping_(stopping),
      method_dt_(),
      method_time_()
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    p | cycle_;
    p | time_;
    p | dt_;
    p | stopping_;
    p | method_dt_;
    p | method_time_;
  };

  //----------------------------------------------------------------------
  /// Initializers
  //----------------------------------------------------------------------

  void set_cycle(int cycle)
  { cycle_ = cycle; }
  void set_time (double time)
  { time_ = time; }
  void set_dt (double dt)
  { dt_ = dt; }
  void set_stopping (bool stopping)
  { stopping_ = stopping; }
  void init (int cycle, double time, double dt, bool stopping)
  {
    set_cycle (cycle);
    set_time (time);
    set_dt (dt);
    set_stopping (stopping);
  }

  void init_method(int n = 0)
  {
    if (n==0) {
      method_dt_.clear();
      method_time_.clear();
    } else {
      method_dt_.resize(n);
      method_time_.resize(n);
      std::fill(method_dt_.begin(),method_dt_.end(), 0.0);
      std::fill(method_time_.begin(),method_time_.end(),0.0);
    }
  }

  //----------------------------------------------------------------------
  /// Accessors
  //----------------------------------------------------------------------

  /// get scalar cycle
  int cycle() const { return cycle_; }
  /// get scalar time
  double time() const { return time_; }
  /// get scalar timestep
  double dt() const { return dt_; }
  /// get scalar stopping criteria
  bool stopping () const { return stopping_; }

  /// get timestep for given method
  double method_dt(int index_method) const {
    ASSERT2("State::method_dt()",
            "dt array length %d is too small for index %d",
            method_dt_.size(),index_method,
            ((0 <= index_method) && (index_method < method_dt_.size())));
    return method_dt_[index_method];
  }
  /// get vector timestep per method
  std::vector<double> & get_method_dt(int n) {
    if (method_dt_.size() < n) method_dt_.resize(n);
    return method_dt_;
  }
  /// get scalar time for given method
  double method_time(int index_method) const {
    ASSERT2("State::method_time()",
            "time array length %d is too small for index %d",
            method_time_.size(),index_method,
            ((0 <= index_method) && (index_method < method_time_.size())));
    return method_time_[index_method];
  }
  /// get vector time per method
  std::vector<double> & get_method_time(int n) {
    if (method_time_.size() < n) method_time_.resize(n);
    return method_time_;
  }

  //----------------------------------------------------------------------
  /// Modifiers
  //----------------------------------------------------------------------

  /// increment cycle by given amount (default one)
  void increment_cycle (int increment = 1)
  { cycle_ += increment; }

  /// increment time by given amount (default dt_)
  void increment_time (double dt = -1)
  { time_ += (dt==-1) ? dt_ : dt; }

private: // functions


protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current global timestep
  double dt_;

  /// Current stopping criteria
  bool stopping_;

  /// Current timestep for each method
  std::vector<double> method_dt_;

  /// Current time for each method
  std::vector<double> method_time_;

};

#endif /* DATA_STATE_HPP */

