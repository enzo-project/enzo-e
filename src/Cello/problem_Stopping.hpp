// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Stopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Method] Declaration of the Stopping class

#ifndef PROBLEM_STOPPING_HPP
#define PROBLEM_STOPPING_HPP

class Stopping : public PUP::able  {

  /// @class    Stopping
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate stopping criteria

public: // interface

  /// Constructor
  Stopping(int    stop_cycle,
	   double stop_time,
	   double stop_seconds) throw()
    : stop_cycle_   (stop_cycle),
      stop_time_    (stop_time),
      stop_seconds_ (stop_seconds)

  {
    DEBUG3 ("cycle %d  time %g  seconds %g\n",  stop_cycle_,stop_time_,stop_seconds_); 
    timer_.start();
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(Stopping);

  /// CHARM++ migration constructor for PUP::able

  Stopping (CkMigrateMessage *m) : PUP::able(m) {}

  /// Destructor
  virtual ~Stopping()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    PUP::able::pup(p);
    // NOTE: change this function whenever attributes change
    p | stop_cycle_;
    p | stop_time_;
    p | stop_seconds_;
  }

  /// Return whether the simulation is done
  virtual bool complete (int    curr_cycle,
			 double curr_time) const throw()
  {
    if ((stop_cycle_ == std::numeric_limits<int>::max()) &&
	(stop_time_  == std::numeric_limits<double>::max()) &&
	(stop_seconds_  == std::numeric_limits<double>::max())) {
      ERROR("Stopping::complete",
	    "No stopping criteria specified");
    }
    bool stop = ( ! ((stop_time_    == -1.0 || curr_time      < stop_time_ ) &&
		     (stop_seconds_ == -1.0 || timer_.value() < stop_seconds_ ) &&
		     (stop_cycle_   == -1   || curr_cycle     < stop_cycle_)));
    
    return stop;
  }

  /// Return stopping cycle
  double stop_cycle () const throw()
  { return stop_cycle_; };

  /// Return stopping time
  double stop_time () const throw()
  { return stop_time_; };

  /// Return stopping seconds
  double stop_seconds () const throw()
  { return stop_seconds_; };

protected:

  /// Stop cycle
  int stop_cycle_;

  /// Stop time
  double stop_time_;

  /// Stop seconds (wall-clock time)
  double stop_seconds_;

  /// Timer
  Timer timer_;

};

#endif /* PROBLEM_STOPPING_HPP */

