// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Method] Declaration of the EnzoStopping class

#ifndef ENZO_ENZO_STOPPING_HPP
#define ENZO_ENZO_STOPPING_HPP

class EnzoStopping : public Stopping {

  /// @class    EnzoStopping
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate stopping criteria

public: // interface

  /// Constructor
  EnzoStopping(int    cycle,
	       double time,
	       double seconds,
	       double redshift) throw()
    : Stopping (cycle, time, seconds),
      stop_redshift_(redshift)

  {
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoStopping);

  /// CHARM++ migration constructor for PUP::able

  EnzoStopping (CkMigrateMessage *m)
    : Stopping (m),
      stop_redshift_(0)
  {}

  /// Destructor
  virtual ~EnzoStopping()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    Stopping::pup(p);
    // NOTE: change this function whenever attributes change
    p | stop_redshift_;
  }

  /// Return whether the simulation is done
  virtual bool complete (int    curr_cycle,
			 double curr_time) const throw()
  {
    bool stop = Stopping::complete(curr_cycle,curr_time);

    // compute current redshift
    double redshift = 0.0;
    return stop || ((stop_redshift_ != -1.0) && (redshift >= stop_redshift_));
  }

  /// Return stopping cycle
  double stop_redshift () const throw()
  { return stop_redshift_; };
protected:

  /// Stop redshift
  double stop_redshift_;

};

#endif /* ENZO_ENZO_STOPPING_HPP */

