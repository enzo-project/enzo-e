// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgState.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2023-04-18
/// @brief    [\ref Charm] Declaration of the MsgState Charm++ Message

#ifndef CHARM_MSG_STATE_HPP
#define CHARM_MSG_STATE_HPP

#include "cello.hpp"

class Simulation;

class MsgState : public CMessage_MsgState {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

//   MsgState();

  virtual ~MsgState();

  MsgState();

  /// Copy constructor
  MsgState(const MsgState & msg_state) throw()
    : CMessage_MsgState() // do NOT call copy constructor on base
  {
    ++counter[cello::index_static()];
    copy_(msg_state);
  }

  MsgState & operator = (const MsgState & msg_state)
  {
    copy_(msg_state);
    return *this;
  }

  /// Copy state variables from this message to the provided Simulation object
  void update (Simulation * simulation);

  void set_state (double time, double dt, int cycle, int stop)
  {
    time_  = time;
    dt_    = dt;
    cycle_ = cycle;
    stop_  = stop;
  }

  void get_state(double & time, double & dt, int & cycle, int & stop)
  {
    time  = time_;
    dt    = dt_;
    cycle = cycle_;
    stop  = stop_;
  }


public: // static methods

  /// Pack data to serialize
  static void * pack (MsgState*);

  /// Unpack data to de-serialize
  static MsgState * unpack(void *);

protected: // methods

  void copy_(const MsgState & msg_state)
  {
    time_  = msg_state.time_;
    dt_    = msg_state.dt_;
    cycle_ = msg_state.cycle_;
    stop_  = msg_state.stop_;
  }

protected: // attributes

  double time_;
  double dt_;
  int cycle_;
  int stop_;

};

#endif /* CHARM_MSG_STATE_HPP */

