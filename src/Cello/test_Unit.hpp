// See LICENSE_CELLO file for license and copyright information

/// @file     test_Unit.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sat Feb 23 15:22:59 PST 2008
/// @brief    [\ref Test] Define the unit namespace and unit testing functions

#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

/// @def      UNIT_MAX_NAME_LEN
/// @brief    Maximum length of a class or function name
#define UNIT_MAX_NAME_LEN 40

/// @def      unit_incomplete
/// @brief    Value used instead of "true" or "false" to indicate an
/// incomplete test
#define unit_incomplete 42

/// @def   unit_assert
/// @brief Assert result of test
#define unit_assert(RESULT)				\
  Unit::instance()->assertion(RESULT, __FILE__,__LINE__)

/// @def   unit_assert_quiet
/// @brief Assert result of test, only output for FAIL
#define unit_assert_quiet(RESULT)				\
  Unit::instance()->assertion(RESULT, __FILE__,__LINE__,true);

/// @def   unit_init
#define unit_init(RANK, SIZE)			\
  Unit::instance()->init(RANK, SIZE)

/// @def   unit_finalize
#define unit_finalize()			\
  Unit::instance()->finalize()

/// @def   unit_set_class
#define unit_class(CLASS_NAME)		\
  Unit::instance()->set_class(CLASS_NAME)

/// @def   unit_set_func
#define unit_func(FUNC_NAME)		\
  Unit::instance()->set_func(FUNC_NAME)

/// @def   unit_set_func
#define unit_func_quiet(CLASS_NAME, FUNC_NAME)		\
  Unit::instance()->set_func(CLASS_NAME,FUNC_NAME)

class Unit {

  /// @class    Unit
  /// @ingroup  Test
  /// @brief    [\ref Test] Class to aid writing unit tests

private:

  /// Private constructor of the Unit object [singleton design pattern]
  Unit();

  /// Private destructor  of the Unit object [singleton design pattern]
  ~Unit();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    WARNING ("Unit::pup","skipping instance_ [static]");

    PUParray (p,class_name_,UNIT_MAX_NAME_LEN);
    PUParray (p,func_name_, UNIT_MAX_NAME_LEN);
    p | test_num_;
    p | is_active_;
    p | comm_size_;
    WARNING ("Unit::pup","pupping comm_rank_");
    p | comm_rank_;
    p | timer_;

  }
#endif


public:

  /// Return an instance of a Unit object
  static Unit * instance()
  { 
    // if ( instance_ == NULL )
    //   instance_ = new Unit;

    // return instance_ ? instance_ : (instance_ = new Unit);
    return & instance_;
  };

  /// Initialize unit testing
  void init (int rank, int size);

  /// Finalize unit testing
  void finalize ();

  /// Set the current unit testing class name
  void set_class (const char * class_name);

  /// Set the current unit testing function name
  void set_func (const char * func_name);

  /// Set the current unit testing class and function names
  void set_func (const char * class_name, const char * func_name);

  /// Assert result of test macro; called by unit_assert macro
  bool assertion (int result, const char * file, int line, bool quiet=false);

private:

  /// Singleton instance of the Unit object
  static Unit instance_;

  /// Output string for passed tests
  static const char * pass_string_;

  /// Output string for failed tests
  static const char * fail_string_;

  /// Output string for incomplete tests
  static const char * incomplete_string_;
  
  /// Name of the current class being tested
  char class_name_ [UNIT_MAX_NAME_LEN];

  /// Name of the current function being tested
  char func_name_  [UNIT_MAX_NAME_LEN];

  /// Running count of tests
  int test_num_;

  /// Whether to output passed tests or not
  int is_active_;

  /// Number of processors
  int comm_size_;

  /// Process rank
  int comm_rank_;

  /// Timer
  Timer timer_;
};

#endif /* TEST_UNIT_HPP */
