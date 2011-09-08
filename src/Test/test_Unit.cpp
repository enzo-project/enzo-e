// See LICENSE_CELLO file for license and copyright information

/// @file     test_Unit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "test.hpp"

//----------------------------------------------------------------------

Unit * Unit::instance_ = 0;

//----------------------------------------------------------------------

Unit::Unit()
  : test_num_(0),
    pass_string_(" pass "),
    fail_string_(" FAIL "),
    incomplete_string_(" incomplete "),
    is_active_(false),
    comm_size_(1),
    comm_rank_(0)
{
}

//----------------------------------------------------------------------

void Unit::init (int rank , int size )
{
  is_active_ = (rank == 0);

  comm_rank_ = rank;
  comm_size_ = size;

  if (is_active_) PARALLEL_PRINTF ("UNIT TEST BEGIN\n");
}

//----------------------------------------------------------------------

void Unit::finalize ()
{
  if (is_active_)  PARALLEL_PRINTF ("UNIT TEST END\n");
}

//----------------------------------------------------------------------

void Unit::set_class (const char * class_name)
{
  strncpy (class_name_, class_name, UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

void Unit::set_func (const char * func_name)
{
  strncpy (func_name_, func_name, UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

void Unit::set_func (const char * class_name, const char * func_name)
{
  strncpy (class_name_,class_name,UNIT_MAX_NAME_LEN);
  strncpy (func_name_, func_name, UNIT_MAX_NAME_LEN);
}

//----------------------------------------------------------------------

bool Unit::assert (int result, const char * file, int line, bool quiet)
{
  // Only print Pass on local process

  if (is_active_ || (result != true) ) { 

    // don't print if quiet is set even if test passed
    if (! (quiet && result)) {
      const char * result_string = 0;
      if (result == true) {
	result_string = pass_string_;
      } else if (result == false) {
	result_string = fail_string_;
      } else {
	result_string = incomplete_string_;
      }
      PARALLEL_PRINTF ("%s %d/%d %s %d %s %s\n",
		       result_string,
		       comm_rank_, 
		       comm_size_,
		       file, 
		       line, 
		       class_name_,
		       func_name_);
      fflush(stdout);
    }
  }
  test_num_++;
  return result;
}
//======================================================================

