// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef TEST_UNIT_HPP
#define TEST_UNIT_HPP

/// @file     test_unit.hpp
/// @brief    Define the unit namespace and unit testing functions
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Convert namespace to class
/// @todo     Address need to edit code to use pass/fail string colors or not
/// @date     Sat Feb 23 15:22:59 PST 2008

/// @def      UNIT_MAX_NAME_LEN
/// @brief    Maximum length of a class or function name
#define UNIT_MAX_NAME_LEN 40

namespace unit {

  /// @namespace unit
  /// @ingroup   Test
  /// @brief     Current class name, function name, and test results
       
  char class_name[UNIT_MAX_NAME_LEN] = {0};
  char func_name [UNIT_MAX_NAME_LEN] = {0};

  int test_num = 1;

  //const char * pass_string = "\033[01;32mPass\033[00m";
  //const char * fail_string = "\033[01;31mFAIL\033[00m";
  const char * pass_string = " [ pass ]";
  const char * fail_string = " [*FAIL*]";

  int process_rank  = 0;
  int process_count = 1;
}

/// @brief Initialize unit testing.  Used to set is_root for parallel runs
void unit_init (int rank = 0, int count = 1)
{
  unit::process_rank  = rank;
  unit::process_count = count;
  unit::class_name[0] = 0;
  unit::func_name[0]  = 0;
  if (unit::process_rank == 0) {
    printf ("UNIT TEST BEGIN\n");
  }
}

/// @brief Finalize unit testing.  Used to ensure test program didn't exit prematurely
void unit_finalize ()
{
  if (unit::process_rank == 0) {
    printf ("UNIT TEST END\n");
  }
}

/// @brief Set the current unit testing class name
void unit_class (const char * c)
{
  strncpy (unit::class_name,c,UNIT_MAX_NAME_LEN);
}

/// @brief Write the size of the class in bytes
template <typename T>
void unit_size ()
{
  if (unit::process_rank == 0) {
    printf ("sizeof (%s) = %lu\n",unit::class_name,sizeof(T));
  }
}

/// @brief Set the current unit testing function name
void unit_func (const char * f)
{
  strncpy (unit::func_name,f,UNIT_MAX_NAME_LEN);
}

/// @brief Assert result of test; macro used for FILE and LINE expansion
#define unit_assert(RESULT) unit_assert_(RESULT, __FILE__,__LINE__);

/// @brief Assert result of test macro; called by unit_assert macro
void unit_assert_ (bool result, const char * file, int line)
{
  if (unit::process_rank == 0 || ! result) {
    printf ("%s %d/%d %s:%d %s::%s() %d\n",
	    (result)? unit::pass_string : unit::fail_string,
	    unit::process_rank, 
	    unit::process_count,
	    file, 
	    line, 
	    unit::class_name,
	    unit::func_name,
	    unit::test_num);
    fflush(stdout);
  }
  unit::test_num++;
}

#endif
