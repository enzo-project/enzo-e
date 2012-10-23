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

  /// Return the number of PAPI events
  int num_events() const throw();

  /// Return the name of the ith event counter
  std::string event_name (int index_event) const throw();

  /// Add a new counter, returning the id
  int add_event(std::string event) throw();

  /// Start event counting for the region
  void start_events() throw();

  /// Stop event counting for the region
  void stop_events() throw();


  /// Return number of regions
  int num_regions() const throw();

  /// Return the current region
  std::string region_name (int index_region) const throw();

  /// Return the index of the given region
  int region_index (std::string name) const throw();

  /// Add a new region, returning the id
  int add_region(std::string region) throw();

  /// Push a new region onto the stack
  void start_region(int index_region) throw();

  /// Push a new region onto the stack
  void stop_region(int index_region) throw();

  /// Clear the counters for the region
  void clear_region(int index_region) throw();

  /// Read the counters for the region
  void read_region(int index_region) throw();

  /// Return array to events for the ith region
  const long long * values (int index_region) const throw();


private: // functions

  void insert_region_(std::string) throw();

private: // attributes


  /// Whether counting has started
  bool is_started_;

  /// Whether PAPI is initialized
  bool is_initialized_;


  /// PAPI event set
  int event_set_;

  /// Number of PAPI events successfully added to the event_set_
  int num_events_;

  /// vector of event names in event set
  std::vector<std::string> event_names_;


  /// Number of regions in lists
  int num_regions_;

  /// list of region names
  std::vector<std::string> region_names_;

  /// list of counter values
  std::vector< std::vector<long long> > region_events_;


  /// mapping of region name to index
  std::map<std::string,int> region_index_;

};

#endif /* PERFORMANCE_PAPI_HPP */

