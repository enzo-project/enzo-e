// See LICENSE_CELLO file for license and copyright information

/// @file     error_Error.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-07
/// @brief    Implementation of Error functions

#include "error.hpp"

//----------------------------------------------------------------------

extern void message2_
(
 FILE *       fp,
 const char * type,
 const char * file,
 int          line,
 const char * function,
 const char * message,
 ...
 )

{
  char buffer[ERROR_LENGTH];

  va_list fargs;
  
  // Process any input arguments

  va_start(fargs,message);
  vsnprintf (buffer,ERROR_LENGTH, message,fargs);
  va_end(fargs);

  Monitor * monitor = Monitor::instance();

  monitor->write (fp,"[%s] %s:%d", type,file,line);
  if (function != "" || buffer != "") {
    monitor->write (fp,"[%s] %s %s",type,function,buffer);
  }
  // monitor->write (fp,"     %10s  %s:%d\n",type,file,line);
  // if (strcmp(function,"") != 0)
  //   monitor->write (fp,"     %10s  %s()\n", type,function);
  // if (strcmp(message,"") != 0)
  //   monitor->write (fp,"     %10s  %s\n",   type,buffer);
  // monitor->write (fp,"\n");

  //  fflush(fp);
}

