// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersUser.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-19
/// @brief    [\ref Lcaperf] Declaration of the CountersUser class

#ifndef LCAPERF_COUNTERS_USER_HPP
#define LCAPERF_COUNTERS_USER_HPP

class CountersUser : public Counters {

  /// @class    CountersUser
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] User counters

public: // interface

  /// Constructor
  CountersUser() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~CountersUser() throw();

  /// Copy constructor
  CountersUser(const CountersUser & counters) throw();

  /// Assignment operator
  CountersUser & operator= (const CountersUser & counters) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    Counters::pup(p);

    p | value_;
  }
#endif

  //----------------------------------------------------------------------

  /// Create a new counter
  void create (std::string counter, counters_type type);

  /// Delete a counter
  void remove (std::string counter);

  /// Increment a user counter
  void increment (std::string counter, long long value);

  /// Assign a value to a user counter
  void assign (std::string counter, long long value);

protected: // functions

  /// Create and start a new set of counters for a new key
  virtual long long * start_();

  /// stop a set of counters
  virtual void stop_(long long * );

  /// Update global counters for the given key
  virtual void update_ (std::string key, long long * counters);
  
  //----------------------------------------------------------------------

private: // functions

  //----------------------------------------------------------------------

private: // attributes

  /// Counter values
  std::vector<long long> value_;

};

#endif /* LCAPERF_COUNTERS_HPP */

