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
  }
#endif

  //----------------------------------------------------------------------
  // global control
  //----------------------------------------------------------------------

  /// Initialize PAPI
  void init() throw();

  /// Start counters
  void start() throw();

  /// Stop counters
  void stop() throw();

  /// Clear the counters
  void clear() throw();

  /// Read the counters
  void read() throw();

  /// Return the name of the ith counter
  std::string name (int id) const throw();

  /// Return the value of the ith counter
  long long value (int id) const throw();

  /// Return the number of counters
  int num_counters() const throw();

  /// Add a new counter, returning the id
  int add_counter(int event) throw();

private: // attributes

  /// Whether counting has started
  bool is_started_;

  /// Whether PAPI is initialized
  bool is_initialized_;

  /// PAPI event set
  int event_set_;

  /// Number of PAPI counters successfully added to the event_set_
  int num_counters_;

  /// list of counter names
  std::vector<std::string> names_;

  /// list of counter values
  std::vector<long long> values_;

  

};

#endif /* PERFORMANCE_PAPI_HPP */

