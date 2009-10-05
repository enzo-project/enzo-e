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
 * @file      TEMPLATE.hpp
 * @brief     Brief description of the file
 * @author    Bart Simpson
 * @date      Thu Feb 21 12:45:56 PST 2008
 * @bug       Example bug description 
 * @note      Modified by Alfred E. Newman, Thu Feb 21 12:49:16 PST 2008
 *            description of modification
 *
 * Detailed description of the file
 *
 * $Id$
 *
 *********************************************************************
 */

#include "performance_timer.hpp"

class Monitor {

/** 
 *********************************************************************
 *
 * @class     Monitor
 * @brief     Brief description of the class
 * @ingroup   Group 
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
  void print (const char * message)
  {
    fprintf (fp_,"%6.1f %s\n",timer_.value(),message);
  };
  
private:

  FILE * fp_;
  Timer timer_; 

};

#endif /* MONITOR_HPP */

