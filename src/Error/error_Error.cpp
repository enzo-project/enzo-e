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

  std::string file_str = file;
  if (strlen(file)>0) {
    // Only print file name not full path
    file_str = file_str.substr(file_str.rfind("/")+1,std::string::npos);
  }

  if (strlen(message) > 0) {
    monitor->write (fp,type,"%s:%d  %s %s",file_str.c_str(),line,function,buffer);
  } else if (strlen(function) > 0) {
    monitor->write (fp,type,"%s:%d  %s",file_str.c_str(),line,function);
  } else {
    monitor->write (fp,type,"%s:%d",file_str.c_str(),line);
  }

  fflush (fp);
}

//----------------------------------------------------------------------

void terminate()
{
  const int buffer_size = 64;

  void * buffer[buffer_size];
  
  int num_symbols = backtrace (buffer,buffer_size);

  char ** symbols = backtrace_symbols(buffer,num_symbols);

  for (int i=0; i<num_symbols; i++) {
    message2_(stderr,"EXIT","",0,"",symbols[i]);
  }
  
  PARALLEL_EXIT;
}
