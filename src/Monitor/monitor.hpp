#ifndef MONITOR_HPP
#define MONITOR_HPP

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

  /// Plot an array as a PNG file
  void plot_png (std::string name, 
		 Scalar * array, int nx, int ny,
		 Scalar min, Scalar max, 
		 const int * color_map, int color_length);
  
private:

  FILE * fp_;
  Timer timer_; 

};

#endif /* MONITOR_HPP */

