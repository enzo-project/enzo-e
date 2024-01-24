// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoState.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2024-01-23
/// @brief    [\ref Enzo] Declaration of the EnzoState class

#ifndef ENZO_STATE_HPP
#define ENZO_STATE_HPP

class EnzoState : public State {

  /// @class    EnzoState
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo]

public: // interface

  /// Constructor
  EnzoState() throw()
  : State(),
    redshift_(0.0)
  {
  }

  EnzoState(int cycle, double time, double dt, bool stopping) throw() :
    State(cycle,time,dt,stopping),
    redshift_(0.0)
  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoState);

  /// CHARM++ migration constructor
  EnzoState(CkMigrateMessage *m)
    : State (m),
      redshift_(0.0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    State::pup(p);

    p | redshift_;
  };

  virtual void set_time (double time)
  {
    time_ = time;
    Simulation * simulation = cello::simulation();
    EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
    EnzoPhysicsCosmology * cosmology =
      units ? units->cosmology() : nullptr ;
    if (cosmology) {
      cosmology->set_current_time(time);
      redshift_ = cosmology->current_redshift();
    }
  }

  void set_redshift (double redshift)
  { redshift_ = redshift; }

  double redshift () const
  { return redshift_; }

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Current redshift
 double redshift_;
};

#endif /* ENZO_STATE_HPP */
