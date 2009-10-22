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


//======================================================================
// LOCAL DEFINES
//======================================================================

#define MEMORY_MAX_NUM_GROUPS 20

//======================================================================
// TYPEDEFS
//======================================================================

typedef int      memory_group_handle;
// typedef unsigned size_t;

class Memory {

/** 
 *********************************************************************
 *
 * @class     Memory
 * @brief     Brief description of the class
 * @ingroup   Group 
 *
 * Detailed description of the class
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

  /// Begin allocating memory associated with the specified group
  static void begin_group ( const char * group_name ) throw ();

  /// End allocating memory associated with the specified group
  static void end_group ( const char * group_name ) throw ();

  /// Return name of the current group
  static const char * current_group () throw ();

  /// Return handle for the current group
  static int current_handle () throw ();

  /// Current number of bytes allocated
  static size_t bytes ( memory_group_handle group_handle = 0 ) throw ();
  
  /// Estimate of amount of local memory availables);
  static size_t available ( memory_group_handle group_handle = 0 ) throw ();

  /// Estimate of used / available memory
  static float efficiency ( memory_group_handle group_handle = 0 ) throw ();

  /// Maximum number of bytes allocated
  static size_t highest ( memory_group_handle group_handle = 0 ) throw ();


  /// Specify the maximum number of bytes to use
  static void set_limit ( size_t size, memory_group_handle group_handle = 0)
    throw ();

  /// Query the maximum number of bytes to use
  static size_t get_limit ( memory_group_handle group_handle = 0 ) throw ();


  /// Return the number of calls to allocate for the group
  static int num_new ( memory_group_handle group_handle = 0 ) throw ();

  /// Return the number of calls to deallocate for the group
  static int num_delete ( memory_group_handle group_handle = 0 ) throw ();

  static void print () throw ();

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  static void check_handle_(memory_group_handle group_handle) 
    throw (ExceptionMemoryBadGroupHandle())
      {  if (!(0 <= group_handle && group_handle < num_groups_)) {
	  throw (ExceptionMemoryBadGroupHandle());
	}
      }
  
//   void *new_(size_t bytes) throw (std::bad_alloc);
//   void delete_(void *p) throw ();

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// The current group index, or 0 if none

  static int curr_group_;

  /// Current number of known groups

  static int num_groups_; 

  /// Array of known group names

  static char * group_names_ [MEMORY_MAX_NUM_GROUPS];

  /// Hardware parameters

  static size_t available_   [MEMORY_MAX_NUM_GROUPS];

  /// Collected statistics for different groups

  static size_t bytes_       [MEMORY_MAX_NUM_GROUPS];
  static size_t bytes_high_  [MEMORY_MAX_NUM_GROUPS];
  static size_t new_calls_   [MEMORY_MAX_NUM_GROUPS];
  static size_t new_bytes_   [MEMORY_MAX_NUM_GROUPS];
  static size_t delete_calls_[MEMORY_MAX_NUM_GROUPS];
  static size_t delete_bytes_[MEMORY_MAX_NUM_GROUPS];
};


#endif /* MEMORY_HPP */
