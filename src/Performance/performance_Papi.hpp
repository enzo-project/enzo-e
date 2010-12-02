// $Id: performance_Papi.hpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PERFORMANCE_PAPI_HPP
#define PERFORMANCE_PAPI_HPP

/// @file     performance_Papi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Dec  2 11:22:30 PST 2010
/// @brief    Interface to the PAPI library


class Papi {

  /// @class    Papi
  /// @ingroup  Groupname
  /// @brief    Brief description of class Papi.

public: // interface

  /// Constructor
  Papi() throw()
    : is_started_(false),
      time_real_total_(0),
      time_proc_total_(0),
      flop_count_total_(0),
      mflop_rate_(0),
      time_real_(0),
      time_proc_(0),
      flop_count_(0)
      
  {};

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Papi() throw()
  {};

  /// Copy constructor
  Papi(const Papi & papi) throw()
  {};

  /// Assignment operator
  Papi & operator= (const Papi & papi) throw()
  {};

  //----------------------------------------------------------------------
  // global control
  //----------------------------------------------------------------------

  /// Initialize counters
  void init() throw()
  {
#ifdef CONFIG_USE_PAPI
    PAPI_flops(&time_real_total_, 
	       &time_proc_total_, 
	       &flop_count_total_,
	       &mflop_rate_);
#endif
  }

  /// Start counters
  void start() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (is_started_) {
      WARNING_MESSAGE("Papi::start",
		      "Counters already started");
    } else {
      is_started_ = true;
      PAPI_flops(&time_real_total_, 
		 &time_proc_total_, 
		 &flop_count_total_,
		 &mflop_rate_);
    }
#endif
  }

  /// Stop counters
  void stop() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (! is_started_) {
      WARNING_MESSAGE("Papi::stop",
		      "Counters already stopped");
    } else {
      is_started_ = false;
      PAPI_flops(&time_real_, 
		 &time_proc_, 
		 &flop_count_,
		 &mflop_rate_);

      time_real_  = time_real_  - time_real_total_;
      time_proc_  = time_proc_  - time_proc_total_;
      flop_count_ = flop_count_ - flop_count_total_;
    }
#endif
  }


  /// Real time between start() and stop()
  float time_real() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (is_started_) {
      WARNING_MESSAGE("Papi::time_real",
		      "Counters must be stopped");
      return 0.0;
    } else {
      return time_real_;
    }
#else
    return 0.0;
#endif
  }

  /// Process time between start() and stop()
  float time_proc() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (is_started_) {
      WARNING_MESSAGE("Papi::time_proc",
		      "Counters must be stopped");
      return 0.0;
    } else {
      return time_proc_;
    }
#else
    return 0.0;
#endif
  }

  /// Return number of flops between start() and stop()
  long long flop_count() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (is_started_) {
      WARNING_MESSAGE("Papi::flop_count",
		      "Counters must be stopped");
      return 0;
    } else {
      return flop_count_;
    }
#else
    return 0;
#endif
  }

  /// Return MFlop rate between start() and stop()
  float mflop_rate() throw()
  {
#ifdef CONFIG_USE_PAPI
    if (is_started_) {
      WARNING_MESSAGE("Papi::mflop_rate",
		      "Counters must be stopped");
      return 0.0;
    } else {
      return mflop_rate_;
    }
#else
    return 0.0;
#endif
  }


private: // attributes

  /// Whether counting has started
  bool is_started_;

  /// Argument 1 to PAPI_flops(): real time
  float time_real_total_;

  /// Argument 2 to PAPI_flops(): process time
  float time_proc_total_;

  /// Argument 3 to PAPI_flops(): floating point operations
  long long flop_count_total_;

  /// Argument 4 to PAPI_flops(): floating point rate in MFlops
  float mflop_rate_;

  /// Real time since last start
  float time_real_;

  /// Process time since last stop
  float time_proc_;

  /// Argument 3 to PAPI_flops(): floating point operations
  long long flop_count_;

};

#endif /* PERFORMANCE_PAPI_HPP */

