/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
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

#ifndef MONITOR_HPP
#define MONITOR_HPP


/** 
 *********************************************************************
 *
 * @file      monitor.hpp
 * @brief     Declaration for Monitor clas
 * @author    James Bordner
 * @date      2009-10-05
 *
 * Declaration for Monitor clas
 *
 * $Id$
 *
 *********************************************************************
 */

#include <string>

#include "cello.h"

#include "performance_timer.hpp"

enum enum_reduce {
  reduce_unknown,
  reduce_min,
  reduce_max,
  reduce_avg,
  reduce_sum
};

class Monitor {

/** 
 *********************************************************************
 *
 * @class     Monitor
 * @brief     Brief description of the class
 * @ingroup   Monitor
 *
 * Detailed description of the class
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// Initialize the Monitor object

  Monitor(FILE * fp = stdout) 
    : fp_(fp)
  {  timer_.start();};

  /// Print a message
  void print (std::string message)
  {
    fprintf (fp_,"%6.1f %s\n",timer_.value(),message.c_str());
  };

  /// Generate a PNG image of an array
  void image (std::string name, 
	      Scalar * array, 
	      int nx, int ny, int nz,  // Array dimensions
	      int nx0, int ny0, int nz0,  // lower indices of subarray
	      int nx1, int ny1, int nz1,  // upper indices of subarray
	      int axis, enum_reduce op_reduce, // Axis along which to project
	      Scalar min, Scalar max,  // Limits for color map
	      const double * color_map, int color_length);  // Color map
  
private:

  FILE * fp_;
  Timer timer_; 

};

#endif /* MONITOR_HPP */

