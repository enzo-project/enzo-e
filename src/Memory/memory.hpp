// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MEMORY_HPP
#define MEMORY_HPP

/// @file     memory.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Sep  3 16:29:56 PDT 2009 
/// @brief    Interface for the Memory class.  Uses the Singleton design pattern.

#include <stack>

/// @def      MEMORY_MAX_NUM_GROUPS
/// @brief    Maximum number of groups for memory allocation tracking
#define MEMORY_MAX_NUM_GROUPS 20

/// @var      typedef int memory_group_handle
/// @brief    Type for opaque group handles
typedef int memory_group_handle;

class Memory {

  /// @class    Memory
  /// @ingroup  Memory
  /// @brief    Maintains memory allocation and deallocation statistics

public: // interface

  /// Get single instance of the Memory object
  static Memory * get_instance() throw ()
  { return & instance_; }

  /// Allocate memory
  void * allocate ( size_t size ) 
    throw (ExceptionMemoryBadAllocate());

  /// De-allocate memory
  void deallocate ( void * pointer ) 
    throw (ExceptionMemoryBadDeallocate());

  /// Assign a name to a group
  void new_group ( memory_group_handle group_id, 
			  const char *        group_name ) throw ();

  /// Begin allocating memory associated with the specified group
  void begin_group ( memory_group_handle group_id ) throw ();

  /// End allocating memory associated with the specified group
  void end_group ( memory_group_handle group_id ) throw ();

  /// Return name of the current group
  const char * current_group () throw ();

  /// Return handle for the current group
  memory_group_handle current_handle () throw ();

  /// Current number of bytes allocated
  long long bytes ( memory_group_handle group_handle = 0 ) throw ();
  
  /// Estimate of amount of local memory available
  long long available ( memory_group_handle group_handle = 0 ) throw ();

  /// Estimate of used / available memory
  float efficiency ( memory_group_handle group_handle = 0 ) throw ();

  /// Maximum number of bytes allocated
  long long highest ( memory_group_handle group_handle = 0 ) throw ();


  /// Specify the maximum number of bytes to use
  void set_limit ( long long           size, 
			  memory_group_handle group_handle = 0) throw ();

  /// Query the maximum number of bytes to use
  long long get_limit ( memory_group_handle group_handle = 0 ) throw ();


  /// Return the number of calls to allocate for the group
  int num_new ( memory_group_handle group_handle = 0 ) throw ();

  /// Return the number of calls to deallocate for the group
  int num_delete ( memory_group_handle group_handle = 0 ) throw ();

  /// Print memory summary
  void print () throw ();

  /// Reset memory counters for the current group
  void reset () throw ();

  /// Set whether memory tracking is active or not
  void set_active (bool is_active) throw ()
#ifdef CONFIG_USE_MEMORY
  { is_active_ = is_active; } ;
#else
  { } ;
#endif

private: // functions

  /// Create the (single) Memory object (singleton design pattern)
  Memory() throw () { initialize_(); };

  /// Delete the (single) Memory object (singleton design pattern)
  ~Memory() throw () {};

  /// Copy the (single) Memory object (singleton design pattern)
  Memory (const Memory &) {};

  /// Assign the (single) Memory object (singleton design pattern)
  Memory & operator = (const Memory & memory) { return *this; };

  /// Initialize the memory component
  void initialize_() throw ();

  ///  Check the group handle, and throw an exception if bad
  void check_handle_(memory_group_handle group_handle) 
    throw (ExceptionMemoryBadGroupHandle())
  {  
#ifdef CONFIG_USE_MEMORY
    if ( 0 > group_handle || group_handle > MEMORY_MAX_NUM_GROUPS) {
      throw (ExceptionMemoryBadGroupHandle());
    }
#endif
  }
  
private: // attributes

  /// Single instance of the Memory object (singleton design pattern)
  static Memory instance_;

#ifdef CONFIG_USE_MEMORY

  /// Whether keeping track of memory statistics is active or not
  bool is_active_;

  /// The current group index, or 0 if none
  std::stack<memory_group_handle> curr_group_;

  /// Array of known group names
  char * group_names_ [MEMORY_MAX_NUM_GROUPS + 1];

  /// Hardware parameters
  long long available_   [MEMORY_MAX_NUM_GROUPS + 1];

  /// Current bytes allocated for different groups
  long long bytes_       [MEMORY_MAX_NUM_GROUPS + 1];

  /// High-water bytes allocated for different groups
  long long bytes_high_  [MEMORY_MAX_NUM_GROUPS + 1];

  /// Number of calls to new for different groups
  long long new_calls_   [MEMORY_MAX_NUM_GROUPS + 1];

  /// Number of calls to delete for different groups
  long long delete_calls_[MEMORY_MAX_NUM_GROUPS + 1];

#endif

};

#endif /* MEMORY_HPP */
