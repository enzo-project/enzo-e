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

  class MethodState {
    friend State;
  public:
    void init() {
      dt_        = 0.0;
      time_      = 0.0;
      num_steps_ = 0;
      step_      = 0;
    }
    void pup (PUP::er &p) {
      p | dt_;
      p | time_;
      p | num_steps_;
      p | step_;
    }

    double dt() const { return dt_; }
    double time() const { return time_; }
    int num_steps() const { return num_steps_; }
    int step() const { return step_; }

    void set_dt(double dt) { dt_ = dt; }
    void set_time(double time) { time_ = time; }
    void set_num_steps(int num_steps) { num_steps_ = num_steps; }
    void set_step(int step) { step_ = step; }

  protected:
    /// Method's timestep
    double dt_;
    /// Method's current time_    
    double time_;
    /// Number of steps expected for method ( > 1 for supercycling)
    int num_steps_;
    /// Number of steps remaining for method
    int step_;
  };

public: // interface

  /// Constructor
  State() throw()
  : cycle_(0),
    time_(0.0),
    dt_(0.0),
    stopping_(false),
    method_state_()
 {
  }

  /// Constructor
  State(int cycle, double time, double dt, bool stopping) throw()
    : cycle_(cycle),
      time_(time),
      dt_(dt),
      stopping_(stopping),
      method_state_()
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
    p | method_state_;
  };

  //----------------------------------------------------------------------
  /// Initializers
  //----------------------------------------------------------------------

  void set_cycle(int cycle) { cycle_ = cycle; }
  void set_time (double time) { time_ = time; }
  void set_dt (double dt) { dt_ = dt; }
  void set_stopping (bool stopping) { stopping_ = stopping; }
  
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
      method_state_.clear();
    } else {
      method_state_.resize(n);
      for (auto & m : method_state_)
        m.init();
    }
  }

  //----------------------------------------------------------------------
  /// Accessors
  //----------------------------------------------------------------------

  int cycle() const { return cycle_; }
  double time() const { return time_; }
  double dt() const { return dt_; }
  bool stopping () const { return stopping_; }

  /// Get ith MethodState
  MethodState & method(int index_method) {
    ASSERT2("State::method()",
            "array length %d is too small for index %d",
            method_state_.size(),index_method,
            ((0 <= index_method) && (index_method < method_state_.size())));
    return method_state_[index_method];
  }

  //----------------------------------------------------------------------
  /// Modifiers
  //----------------------------------------------------------------------

  /// increment cycle by given amount (default one)
  void increment_cycle (int increment = 1)
  { cycle_ += increment; }

  /// increment time by given amount (default dt_)
  void increment_time (double dt)
  { time_ += dt; }

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

  /// Method-specific state scalars
  std::vector<MethodState> method_state_;
};

#endif /* DATA_STATE_HPP */

