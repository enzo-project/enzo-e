// $Id: lcaperf_counters_papi.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE file for license and copyright information

#ifndef LCAPERF_COUNTERS_PAPI_HPP
#define LCAPERF_COUNTERS_PAPI_HPP

/// @file     lcaperf_counters_papi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon May 23 11:12:54 PDT 2011
/// @brief    [\ref lcaperf] Declaration of the CountersPapi class

/// @enum
/// @brief Indices of papi counters

class CountersPapi : public CountersUser {

  /// @class    CountersPapi
  /// @ingroup  lcaperf
  /// @brief    [\ref lcaperf] Papi counters

public: // interface

  /// Constructor
  CountersPapi() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~CountersPapi() throw();

  /// Copy constructor
  CountersPapi(const CountersPapi & counters) throw();

  /// Assignment operator
  CountersPapi & operator= (const CountersPapi & counters) throw();

  //----------------------------------------------------------------------

protected: // functions

  /// Create and start a new set of counters for the current key
  virtual long long * start_();

  /// stop a set of counters for the current key
  virtual void stop_(long long * );

  /// Update global counters for the given key
  virtual void update_ (std::string key, long long * counters);
  
  //----------------------------------------------------------------------

private: // functions

  void papi_clear_eventset_();
  void papi_create_eventset_();
  void papi_delete_eventset_();
  int papi_insert_counter_ (char * counter_name);
  int papi_delete_counter_(char * counter_name);
  void papi_start_counters_ ();
  void papi_read_counters_ ();
  void papi_stop_counters_ ();

private: // attributes

  int is_papi_active_;

  int event_set_;                               // handle to PAPI event set
  long long * papi_counters_;
  long long vtime_begin_;

};

#endif /* LCAPERF_COUNTERS_PAPI_HPP */

