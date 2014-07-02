// See LICENSE_CELLO file for license and copyright information

/// @file     memory_Memory.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Sep  3 16:29:56 PDT 2009 
/// @brief    [\ref Memory] Interface for the Memory class.  Uses the
/// Singleton design pattern.

#ifndef MEMORY_MEMORY_HPP
#define MEMORY_MEMORY_HPP

/// @def      MEMORY_MAX_NUM_GROUPS
/// @brief    Maximum number of groups for memory allocation tracking
#define MEMORY_MAX_NUM_GROUPS 20

/// @var      typedef int memory_group_handle
/// @brief    Type for opaque group handles
typedef int memory_group_handle;

class Memory {

  /// @class    Memory
  /// @ingroup  Memory
  /// @brief    [\ref Memory] Manage memory allocation and deallocation

public: // interface

/// Get single instance of the Memory object
  static Memory * instance() throw ()
  {
#ifdef CONFIG_USE_MEMORY
    return & instance_; 
#else
    return 0;
#endif
  }

private: // interface
  /// Create the (single) Memory object (singleton design pattern)
  Memory() throw ()
#ifdef CONFIG_USE_MEMORY
 : is_active_(false)
#endif
  { initialize_(); };

  /// Copy the (single) Memory object (singleton design pattern)
  Memory (const Memory &);

  /// Assign the (single) Memory object (singleton design pattern)
  Memory & operator = (const Memory & memory);

  /// Delete the Memory object
  ~Memory() throw ();

public: // interface

  void pup(PUP::er &p)
  {
    p | is_active_;
    p | do_allocate_fill_; 
    p | allocate_fill_value_;
    p | do_deallocate_fill_; 
    p | deallocate_fill_value_;
    p | max_group_id_;
    WARNING ("Memory::pup()","Skipping curr_group_");
    //    p | curr_group_;
    //    WARNING ("Memory::pup()","Skipping group_names_");
    PUParray(p,group_names_ ,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,limit_,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,bytes_,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,bytes_high_,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,bytes_highest_,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,new_calls_,MEMORY_MAX_NUM_GROUPS + 1);
    PUParray(p,delete_calls_,MEMORY_MAX_NUM_GROUPS + 1);
  }

  /// Allocate memory
  void * allocate ( size_t size ) throw ();

  /// De-allocate memory
  void deallocate ( void * pointer ) throw ();

  /// Assign a name to a group
  void new_group ( memory_group_handle group_id, 
		   const char *        group_name ) throw ();

  /// Begin allocating memory associated with the specified group
  void begin_group ( memory_group_handle group_id ) throw ();

  /// End allocating memory associated with the specified group
  void end_group ( memory_group_handle group_id ) throw ();

  /// Return handle and name of the current group
  memory_group_handle current_group (const char ** ) throw ();

  /// Current number of bytes allocated
  long long bytes ( memory_group_handle group_handle = 0 ) throw ();
  
  /// Estimate of used / available memory
  float efficiency ( memory_group_handle group_handle = 0 ) throw ();

  /// Maximum number of bytes allocated within an interval
  long long bytes_high ( memory_group_handle group_handle = 0 ) throw ();

  /// Maximum number of bytes allocated during run
  long long bytes_highest ( memory_group_handle group_handle = 0 ) throw ();

  /// Specify the maximum number of bytes to use
  void set_limit ( long long           size, 
		   memory_group_handle group_handle = 0) throw ();

  /// Return the maximum number of bytes to use
  long long limit ( memory_group_handle group_handle = 0) throw ();

  /// Query the maximum number of bytes left to use, 0 if unknown
  long long available ( memory_group_handle group_handle = 0 ) throw ();


  /// Return the number of calls to allocate for the group
  int num_new ( memory_group_handle group_handle = 0 ) throw ();

  /// Return the number of calls to deallocate for the group
  int num_delete ( memory_group_handle group_handle = 0 ) throw ();

  /// Print memory summary
  void print () throw ();

  /// Reset memory counters for the current group
  void reset () throw ();

  /// Reset bytes_high to current
  void reset_high () throw ();

  /// Set whether memory tracking is active or not
  void set_active (bool is_active) throw ()
#ifdef CONFIG_USE_MEMORY
  { is_active_ = is_active; } ;
#else
  { } ;
#endif

  bool is_active () const throw()
  { 
#ifdef CONFIG_USE_MEMORY
    return is_active_; 
#else
    return 0;
#endif
  }

  /// Set whether to fill memory with a value after allocating
  void set_allocate_fill (bool do_fill, char fill_value = 0)
  {
#ifdef CONFIG_USE_MEMORY
    do_allocate_fill_    = do_fill;
    allocate_fill_value_ = fill_value;
#endif
  };

  /// Set whether to fill memory with a value before deallocating
  void set_deallocate_fill (bool do_fill, char fill_value = 0)
  {
#ifdef CONFIG_USE_MEMORY
    do_deallocate_fill_    = do_fill;
    deallocate_fill_value_ = fill_value;
#endif
  };

private: // functions

  /// Initialize the memory component
  void initialize_() throw ();

  /// Finalize the memory component
  void finalize_() throw ();

  ///  Check the group handle
  void check_handle_(memory_group_handle group_handle) throw ();

private: // attributes

#ifdef CONFIG_USE_MEMORY
  /// Single instance of the Memory object (singleton design pattern)
  static Memory instance_;
#endif

#ifdef CONFIG_USE_MEMORY
  
  /// Whether keeping track of memory statistics is active or not
  bool is_active_;

  /// Whether to fill memory after allocation
  bool do_allocate_fill_; 

  /// allocate clear value
  char allocate_fill_value_;

  /// Whether to clear memory before deallocate
  bool do_deallocate_fill_; 

  /// deallocate clear value
  char deallocate_fill_value_;

  /// The highest active group index
  int max_group_id_;

  /// The current group index, or 0 if none
  std::stack<memory_group_handle> curr_group_;

  /// Array of known group names
  std::string group_names_ [MEMORY_MAX_NUM_GROUPS + 1];

  /// Limit on number of bytes to allocate.  Currently not checked.
  long long limit_   [MEMORY_MAX_NUM_GROUPS + 1];

  /// Current bytes allocated for different groups
  long long bytes_       [MEMORY_MAX_NUM_GROUPS + 1];

  /// Intervaled high-water bytes allocated for different groups
  long long bytes_high_  [MEMORY_MAX_NUM_GROUPS + 1];

  /// High-water bytes allocated for different groups
  long long bytes_highest_  [MEMORY_MAX_NUM_GROUPS + 1];

  /// Number of calls to new for different groups
  long long new_calls_   [MEMORY_MAX_NUM_GROUPS + 1];

  /// Number of calls to delete for different groups
  long long delete_calls_[MEMORY_MAX_NUM_GROUPS + 1];

#endif

};

#endif /* MEMORY_MEMORY_HPP */
