// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Stopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Method] Declaration of the Stopping class

#ifndef PROBLEM_STOPPING_HPP
#define PROBLEM_STOPPING_HPP

class Stopping {

  /// @class    Stopping
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate stopping criteria

public: // interface

  /// Constructor
  Stopping(int    stop_cycle = std::numeric_limits<int>::max(),
	   double stop_time  = std::numeric_limits<double>::max()) throw()
    : stop_cycle_(stop_cycle),
      stop_time_ (stop_time)
  {    DEBUG2 ("cycle %d    time  %g\n",  stop_cycle_,stop_time_); }

  /// Destructor
  virtual ~Stopping()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | stop_cycle_;
    p | stop_time_;

  }

  /// Return whether the simulation is done
  virtual bool complete (int    curr_cycle,
			 double curr_time) const throw()
  {
    if ((stop_cycle_ == std::numeric_limits<int>::max()) &&
	(stop_time_  == std::numeric_limits<double>::max())) {
      ERROR("Stopping::complete",
	    "Neither Stopping::time_stop nor Stopping::cycle_stop initialized");
    }
    bool stop = ( ! ((stop_time_  == -1.0 || curr_time  < stop_time_ ) &&
		     (stop_cycle_ == -1   || curr_cycle < stop_cycle_)));
    //    printf ("DEBUG cycle %d %d   time %g %g  stop %d\n",
    //	    curr_cycle,stop_cycle_,curr_time,stop_time_,stop);

    
    return stop;
  }

  /// Return stopping cycle
  double stop_cycle () const throw()
  { return stop_cycle_; };

  /// Return stopping time
  double stop_time () const throw()
  { return stop_time_; };

protected:

  /// Stop cycle
  int stop_cycle_;

  /// Stop time
  double stop_time_;

};

#endif /* PROBLEM_STOPPING_HPP */

