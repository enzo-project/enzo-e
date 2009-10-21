#ifndef COUNTERS_HPP
#define COUNTERS_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      counters.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include "memory.hpp"

class Counters {

/** 
 *********************************************************************
 *
 * @class     Counters
 * @brief     
 * @ingroup   Performance
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// 
  Counters(size_t num_attributes, size_t num_counters)
    {
      Memory::begin_group("Performance");
      a_  = new int       [num_attributes];
      c_  = new long long [num_counters];
      dc_ = new long long [num_counters];
      Memory::end_group("Performance");
    }

  /// 
  ~Counters()
    {
      Memory::begin_group("Performance");
      delete [] a_;
      delete [] c_;
      delete [] dc_;
      Memory::end_group("Performance");
    }


private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Array of attribute values
  int       * a_;

  /// Array of the current counter values at the start of the region
  long long * c_;

  /// Array of changes to the counter values
  long long * dc_;

};

#endif /* COUNTERS_HPP */

