// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_CountersBasic.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-05-19
/// @brief    [\ref Lcaperf] Declaration of the CountersBasic class
///

#ifndef LCAPERF_COUNTERS_BASIC_HPP
#define LCAPERF_COUNTERS_BASIC_HPP

namespace lca {

class CountersBasic : public Counters {

  /// @class    CountersBasic
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Basic counters

public: // interface

  /// Constructor
  CountersBasic() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~CountersBasic() throw();

  /// Copy constructor
  CountersBasic(const CountersBasic & CountersBasic) throw();

  /// Assignment operator
  CountersBasic & operator= (const CountersBasic & CountersBasic) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    Counters::pup(p);
  }
#endif

  //----------------------------------------------------------------------

protected: // functions

  /// Create and start a new set of counters for a new key
  virtual long long * start_();

  /// stop a set of counters
  virtual void stop_(long long * );

  /// Update global counters for the given key
  virtual void update_ (std::string key, long long * counters);
  
  //----------------------------------------------------------------------

private: // functions

  mutable long long time_begin_;

  inline long long wall_time_() const
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday (&tv,&tz);
    long long value = 1000000L * tv.tv_sec + tv.tv_usec;
    if (time_begin_==-1) time_begin_=value;
    return value - time_begin_;
  };

};

}

#endif /* LCAPERF_COUNTERS_BASIC_HPP */

