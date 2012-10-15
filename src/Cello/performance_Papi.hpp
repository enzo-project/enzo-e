// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Papi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Dec  2 11:22:30 PST 2010
/// @brief    [\ref Performance] Interface to the PAPI library

#ifndef PERFORMANCE_PAPI_HPP
#define PERFORMANCE_PAPI_HPP

class Papi {

  /// @class    Papi
  /// @ingroup  Groupname
  /// @brief    [\ref Performance] Class for accessing PAPI counters

public: // interface

  /// Constructor
  Papi() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  // ~Papi() throw()
  // {};

  // /// Copy constructor
  // Papi(const Papi & papi) throw()
  // {};

  // /// Assignment operator
  // Papi & operator= (const Papi & papi) throw()
  // {return *this};

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    WARNING("PerformancePapi::pup","skipping");
    return;
    p | is_started_;
    p | time_real_total_;
    p | time_proc_total_;
    p | flop_count_total_;
    p | flop_rate_;
    p | time_real_;
    p | time_proc_;
    p | flop_count_;
  }
#endif

  //----------------------------------------------------------------------
  // global control
  //----------------------------------------------------------------------

  /// Start counters
  void start() throw();

  /// Stop counters
  void stop() throw();

  /// Update counters
  void update() throw()
  {  stop(); start(); }

  /// Real time between start() and stop()
  float time_real() const throw();

  /// Process time between start() and stop()
  float time_proc() const throw();

  /// Return number of flops between start() and stop()
  long long flop_count() const throw();

  /// Return flop rate between start() and stop()
  float flop_rate() const throw();

  // void print () const throw();

private: // attributes

  /// Whether counting has started
  bool is_started_;

  /// Argument 1 to PAPI_flops(): real time
  float time_real_total_;

  /// Argument 2 to PAPI_flops(): process time
  float time_proc_total_;

  /// Argument 3 to PAPI_flops(): floating point operations
  long long flop_count_total_;

  /// Argument 4 to PAPI_flops(): floating point rate in flops/s
  float flop_rate_;

  /// Real time since last start
  float time_real_;

  /// Process time since last stop
  float time_proc_;

  /// Argument 3 to PAPI_flops(): floating point operations
  long long flop_count_;

};

#endif /* PERFORMANCE_PAPI_HPP */

