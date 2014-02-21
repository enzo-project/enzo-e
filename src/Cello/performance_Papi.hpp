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
  /// @brief    [\ref Performance] Class for accessing PAPI events

public: // interface

  /// Constructor
  Papi() throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    WARNING("PerformancePapi::pup","skipping");
    return;
  }

  //----------------------------------------------------------------------
  // global control
  //----------------------------------------------------------------------

  /// Initialize PAPI
  void init() throw();

  /// Return the number of PAPI events
  int num_events() const throw();

  /// Return the name of the ith event counter
  std::string event_name (int index_event) const throw();

  /// Add a new counter, returning the id
  int add_event(std::string event) throw();

  /// Start event counting
  void start_events() throw();

  /// Stop event counting
  void stop_events() throw();

  /// Return array to events
  int event_values (long long * values) const throw();


private: // attributes

  /// Whether PAPI is initialized
  bool is_initialized_;

  /// Whether counting has started
  bool is_started_;


  /// PAPI event set
  int event_set_;

  /// Number of PAPI events successfully added to the event_set_
  int num_events_;

  /// vector of event names in event set
  std::vector<std::string> event_names_;

};

#endif /* PERFORMANCE_PAPI_HPP */

