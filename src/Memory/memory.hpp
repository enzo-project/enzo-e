#ifndef MEMORY_HPP
#define MEMORY_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      memory.hpp
 * @brief     Functions for dynamic memory management
 * @author    James Bordner
 * @date      Thu Sep  3 16:29:56 PDT 2009
 * @bug       
 * @note      
 *
 * Functions for dynamic memory management
 *
 * $Id$
 *
 *********************************************************************
 */

typedef unsigned memory_group_handle;

#ifdef USE_MEMORY

//======================================================================
// LOCAL DEFINES
//======================================================================

#define MEMORY_MAX_NUM_GROUPS 20

//======================================================================
// TYPEDEFS
//======================================================================

class Memory {

/** 
 *********************************************************************
 *
 * @class     Memory
 * @brief     Maintains memory allocation and deallocation statistics
 * @ingroup   Group 
 *
 * Maintains memory allocation and deallocation statistics
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// Initialize the memory component
  Memory() throw ();

  /// Allocate memory
  static void * allocate ( size_t size ) 
    throw (ExceptionMemoryBadAllocate());

  /// De-allocate memory
  static void deallocate ( void * pointer ) 
    throw (ExceptionMemoryBadDeallocate());

  /// Assign a name to a group
  static void new_group ( memory_group_handle group_id, const char * group_name ) throw ();

  /// Begin allocating memory associated with the specified group
  static void begin_group ( memory_group_handle group_id ) throw ();

  /// End allocating memory associated with the specified group
  static void end_group ( memory_group_handle group_id ) throw ();

  /// Return name of the current group
  static const char * current_group () throw ();

  /// Return handle for the current group
  static memory_group_handle current_handle () throw ();

  /// Current number of bytes allocated
  static long long bytes ( memory_group_handle group_handle = 0 ) throw ();
  
  /// Estimate of amount of local memory availables);
  static long long available ( memory_group_handle group_handle = 0 ) throw ();

  /// Estimate of used / available memory
  static float efficiency ( memory_group_handle group_handle = 0 ) throw ();

  /// Maximum number of bytes allocated
  static long long highest ( memory_group_handle group_handle = 0 ) throw ();


  /// Specify the maximum number of bytes to use
  static void set_limit ( long long size, memory_group_handle group_handle = 0)
    throw ();

  /// Query the maximum number of bytes to use
  static long long get_limit ( memory_group_handle group_handle = 0 ) throw ();


  /// Return the number of calls to allocate for the group
  static int num_new ( memory_group_handle group_handle = 0 ) throw ();

  /// Return the number of calls to deallocate for the group
  static int num_delete ( memory_group_handle group_handle = 0 ) throw ();

  static void print () throw ();

  static void reset () throw ();

  static void set_active (bool is_active) throw ()
  { is_active_ = is_active; } ;

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  static void check_handle_(memory_group_handle group_handle) 
    throw (ExceptionMemoryBadGroupHandle())
      {  if ( group_handle > MEMORY_MAX_NUM_GROUPS) {
	  throw (ExceptionMemoryBadGroupHandle());
	}
      }
  
//   void *new_(size_t bytes) throw (std::bad_alloc);
//   void delete_(void *p) throw ();

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------


  /// Whether keeping track of memory statistics is active or not
  static bool is_active_;

  /// The current group index, or 0 if none

  static memory_group_handle curr_group_;

  /// Array of known group names

  static char * group_names_ [MEMORY_MAX_NUM_GROUPS + 1];

  /// Hardware parameters

  static long long available_   [MEMORY_MAX_NUM_GROUPS + 1];

  /// Collected statistics for different groups

  static long long bytes_       [MEMORY_MAX_NUM_GROUPS + 1];
  static long long bytes_high_  [MEMORY_MAX_NUM_GROUPS + 1];
  static long long new_calls_   [MEMORY_MAX_NUM_GROUPS + 1];
  static long long delete_calls_[MEMORY_MAX_NUM_GROUPS + 1];
};

#else /* USE_MEMORY */

class Memory {

/** 
*********************************************************************
*
* @class     Memory
* @brief     Dummy functions for not maintaining memory statistics
* @ingroup   Group 
*
* Dummy functions for not maintaining memory statistics
*
*********************************************************************
*/

public:

//-------------------------------------------------------------------
// PUBLIC OPERATIONS
//-------------------------------------------------------------------

/// Initialize the memory component
Memory() throw () {};

/// Allocate memory
static void * allocate ( size_t size ) { return 0; };

/// De-allocate memory
static void deallocate ( void * pointer ) {};

/// Assign a name to a group
static void new_group ( memory_group_handle, const char *) {};

/// Begin allocating memory associated with the specified group
static void begin_group ( memory_group_handle group_id ) {};

/// End allocating memory associated with the specified group
static void end_group ( memory_group_handle group_id ) {};

/// Return name of the current group
static const char * current_group () {return ""; };

/// Return handle for the current group
static memory_group_handle current_handle () {return 0; };

/// Current number of bytes allocated
static long long bytes ( memory_group_handle group_handle = 0 ) {return 0; };
  
/// Estimate of amount of local memory availables);
static long long available ( memory_group_handle group_handle = 0 ) {return 0; };

/// Estimate of used / available memory
static float efficiency ( memory_group_handle group_handle = 0 ) {return 0; };

/// Maximum number of bytes allocated
static long long highest ( memory_group_handle group_handle = 0 ) {return 0; };


/// Specify the maximum number of bytes to use
static void set_limit ( long long size, memory_group_handle group_handle = 0)
  {};

/// Query the maximum number of bytes to use
static long long get_limit ( memory_group_handle group_handle = 0 ) {return 0; };


/// Return the number of calls to allocate for the group
static int num_new ( memory_group_handle group_handle = 0 ) {return 0; };

/// Return the number of calls to deallocate for the group
static int num_delete ( memory_group_handle group_handle = 0 ) {return 0; };

static void print () {};

static void reset () {};

static void set_active (bool is_active) {}
};

#endif

#endif /* MEMORY_HPP */
