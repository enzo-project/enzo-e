// $Id: lcaperf_counters_mem.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE file for license and copyright information

#ifndef LCAPERF_COUNTERS_MEM_HPP
#define LCAPERF_COUNTERS_MEM_HPP

/// @file     lcaperf_counters_mem.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon May 23 11:12:54 PDT 2011
/// @brief    [\ref lcaperf] Declaration of the CountersMem class

/// @enum
/// @brief Indices of memory counters

class CountersMem : public CountersUser {

  /// @class    CountersMem
  /// @ingroup  lcaperf
  /// @brief    [\ref lcaperf] Mem counters

public: // interface

  /// Constructor
  CountersMem() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~CountersMem() throw();

  /// Copy constructor
  CountersMem(const CountersMem & counters) throw();

  /// Assignment operator
  CountersMem & operator= (const CountersMem & counters) throw();

  //----------------------------------------------------------------------

public: // functions

  static void * my_new (size_t bytes) throw (std::bad_alloc);
  static void   my_delete (void *p) throw ();

  //----------------------------------------------------------------------

protected: // functions

  /// Create and start a new set of counters for the current key
  virtual long long * start_();

  /// stop a set of counters for the current key
  virtual void stop_(long long * );

  /// Update global counters for the given key
  // virtual void update_ (std::string key, long long * counters);

  /// Function for copying static mem counters to lcaperf user counters
  void assign_all_();

private: // attributes

  static long long curr_bytes_;
  static long long high_bytes_;
  static long long new_count_;
  static long long new_bytes_;
  static long long del_count_;
  static long long del_bytes_;

};

#endif /* LCAPERF_COUNTERS_MEM_HPP */

