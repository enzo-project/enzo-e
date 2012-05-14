// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_Counters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-19
/// @brief    [\ref Lcaperf] Declaration of the Counters class

#ifndef LCAPERF_COUNTERS_HPP
#define LCAPERF_COUNTERS_HPP

namespace lca {

//----------------------------------------------------------------------

enum counter_type {
  counter_type_unknown  = 0,
  counter_type_relative = 1,
  counter_type_absolute = 2
};

//----------------------------------------------------------------------

class Counters {

  /// @class    Counters
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Base class for various performance counters

  friend class ItCounterKeys;
  friend class LcaPerf;

public: // interface

  /// Constructor
  Counters(int num_counters) throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Counters() throw();

  /// Copy constructor
  Counters(const Counters & Counters) throw();

  /// Assignment operator
  Counters & operator= (const Counters & Counters) throw();

  /// Initialize the Counters object
  virtual void initialize()
  {};

  /// Begin logging counters to disk
  void begin ();

  /// End logging counters to disk
  void end ();

  /// Start counters for a region
  void start (std::string region,
	      const Attributes * attributes);

  /// Stop counters for a region
  void stop (std::string region,
	     const Attributes * attributes);

  /// Return the value for a given counter
  long long value (std::string key, std::string counter);

  /// Clear all counter values
  void clear ();

  /// Display summary of current counter values
  void print (FILE * fp = stdout);

  /// Return index of the given counter
  int index (std::string counter) const
  { return index_.at(counter); }

  /// Return name for the given counter index
  std::string name (int index) const
  { return name_.at(index); }

  /// Return the user counter type
  counter_type type (int index) const 
  { return type_.at(index); }
  

  /// Return the number of counters
  int num_counters () const 
  { return num_counters_; };

  //----------------------------------------------------------------------

protected: // virtual functions

  /// Generate a key given region and current attribute values
  std::string generate_key_ (std::string        region, 
			     const Attributes * attributes) const;

  /// Merge keys for attributes that may have changed
  std::string merge_keys_ (std::string key1, std::string key2) const;

  //----------------------------------------------------------------------

protected: // virtual functions

  /// Create and start a new set of counters for a new key
  virtual long long * start_() = 0;

  /// stop a set of counters
  virtual void stop_(long long * counters) = 0;

  /// Update global counters for the given key
  virtual void update_ (std::string key, long long * counters) = 0;
  
  //----------------------------------------------------------------------

protected: // attributes

  /// Number of counters
  int num_counters_;

  /// Whether the counters are active
  bool is_tracing_active_;

  /// Whether the counters are active
  bool is_logging_active_;

  /// Global counter values for all keys
  std::map<std::string,long long *> global_;

  /// Map counter name to counter index
  std::map<std::string,int> index_;

  /// Map counter index to counter name
  std::vector<std::string> name_;

  /// Couner types: i.e. absolute or relative
  std::vector<enum counter_type> type_;

  /// Frame stack
  std::stack<std::string> frame_;

  /// Counters stack
  std::stack<long long *> counters_;
};

}

#endif /* LCAPERF_COUNTERS_HPP */

