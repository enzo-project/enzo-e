// See LICENSE_CELLO file for license and copyright information

/// @file     test_Unit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "test.hpp"

//----------------------------------------------------------------------

// Unit * Unit::instance_ = 0;
Unit Unit::instance_;
const char * Unit::pass_string_       = " pass ";
const char * Unit::fail_string_       = " FAIL ";
const char * Unit::incomplete_string_ = " incomplete ";

//----------------------------------------------------------------------

Unit::Unit()
  : test_num_(0),
    is_active_(false),
    comm_size_(1),
    comm_rank_(0)
{
}

//----------------------------------------------------------------------

Unit::~Unit()
{
}

//----------------------------------------------------------------------

void Unit::init (int rank , int size )
{
  is_active_ = (rank == 0);

  timer_.start();

  comm_rank_ = rank;
  comm_size_ = size;

  if (is_active_) {
    PARALLEL_PRINTF ("UNIT TEST BEGIN\n");
    fflush(stdout);
  }
}

//----------------------------------------------------------------------

void Unit::finalize ()
{
  timer_.stop();
  if (is_active_)  {
    PARALLEL_PRINTF ("UNIT TEST END %f\n",timer_.value());
    fflush(stdout);
  }
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

bool Unit::assertion (int result, const char * file, int line, bool quiet)
{
  // Only print Pass on local process

  if (is_active_ || (result != true) ) { 

    // don't print if quiet is set even if test passed
    if (! (quiet && result)) {
      const char * result_string = 0;
      if (result == true) {
	result_string = Unit::pass_string_;
      } else if (result == false) {
	result_string = Unit::fail_string_;
      } else {
	result_string = Unit::incomplete_string_;
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

