//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef COUNTERS_HPP
#define COUNTERS_HPP

///
/// @brief     
/// @author    
/// @date      
/// @ingroup
/// @note      
///

/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */


/** 
 *********************************************************************
 *
 * @file      counters.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello.h"

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
      Memory::begin_group(component_performance);
      a_  = new int       [num_attributes];
      c_  = new long long [num_counters];
      dc_ = new long long [num_counters];
      Memory::end_group(component_performance);
    }

  /// 
  ~Counters()
    {
      Memory::begin_group(component_performance);
      delete [] a_;
      delete [] c_;
      delete [] dc_;
      Memory::end_group(component_performance);
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

