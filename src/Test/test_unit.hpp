// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

/// @file     test_unit.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Convert namespace to class
/// @todo     Address need to edit code to use pass/fail string colors or not
/// @todo     Use all uppercase for visibility, e.g. UNIT_ASSERT(...)
/// @date     Sat Feb 23 15:22:59 PST 2008
/// @brief    [\ref Test] Define the unit namespace and unit testing functions

/// @def      UNIT_MAX_NAME_LEN
/// @brief    Maximum length of a class or function name
#define UNIT_MAX_NAME_LEN 40

#define unit_incomplete 42
namespace unit {

  /// @namespace unit
  /// @ingroup   Test
  /// @brief     [\ref Test] Current class name, function name, and test results
       
  char class_name [UNIT_MAX_NAME_LEN] = {"CLASS_NAME"};
  char func_name  [UNIT_MAX_NAME_LEN] = {"FUNC_NAME"};

  int test_num = 0;

  //const char * pass_string = "\033[01;32mPass\033[00m";
  //const char * fail_string = "\033[01;31mFAIL\033[00m";
  const char * pass_string = " pass ";
  const char * fail_string = " FAIL ";
  const char * incomplete_string = " incomplete ";

  int process_rank  = 0;
  int process_count = 1;
}

/// @brief Initialize unit testing.  Used to set is_root for parallel runs
void unit_init (int rank = 0, int count = 1)
{
  unit::process_rank  = rank;
  unit::process_count = count;
  if (unit::process_rank == 0) {
    PARALLEL_PRINTF ("UNIT TEST BEGIN\n");
  }
}

/// @brief Finalize unit testing.  Used to ensure test program didn't exit prematurely
void unit_finalize ()
{
  if (unit::process_rank == 0) {
    PARALLEL_PRINTF ("UNIT TEST END\n");
  }
}

/// @brief Set the current unit testing function name
void unit_func (const char * class_name, const char * func_name)
{
  strncpy (unit::class_name,class_name,UNIT_MAX_NAME_LEN);
  strncpy (unit::func_name, func_name, UNIT_MAX_NAME_LEN);
}

/// @brief Write the size of the class in bytes
template <typename T>
void unit_size ()
{
  if (unit::process_rank == 0) {
    PARALLEL_PRINTF ("sizeof (%s) = %lu\n",unit::class_name,sizeof(T));
  }
}

/// @brief Assert result of test; macro used for FILE and LINE expansion
#define unit_assert(RESULT) \
  unit_assert_(RESULT, __FILE__,__LINE__);

/// @brief Assert result of test, only output for FAIL; macro used for FILE and LINE expansion
#define unit_assert_quiet(RESULT) \
  unit_assert_(RESULT, __FILE__,__LINE__,true);

/// @brief Assert result of test macro; called by unit_assert macro
bool unit_assert_ (int result, const char * file, int line, bool quiet=false)
{
  // Only print Pass on local process
  if (unit::process_rank == 0 || ! result) { 
    // don't print if quiet is set and test passed
    if (! ( quiet & result)) {
      const char * result_string = 0;
      if (result == true) {
	result_string = unit::pass_string;
      } else if (result == false) {
	result_string = unit::fail_string;
      } else {
	result_string = unit::incomplete_string;
      }
      PARALLEL_PRINTF ("%s %d/%d %s %d %s %s\n",
		       result_string,
		       unit::process_rank, 
		       unit::process_count,
		       file, 
		       line, 
		       unit::func_name,
		       unit::func_name);
      fflush(stdout);
    }
  }
  unit::test_num++;
  return result;
}

#endif
