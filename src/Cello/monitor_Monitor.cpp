// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"

#include "monitor.hpp" 

#include "auto_config.def"


//----------------------------------------------------------------------
Monitor Monitor::instance_; // singleton design pattern)
//----------------------------------------------------------------------

Monitor::Monitor()
  : timer_(new Timer),
    active_(true)
{ 
  timer_->start();

  // Default: always output
  group_default_ = true;

  // turn off debugging
  group_active_["DEBUG"] = false;
}

//----------------------------------------------------------------------

Monitor::~Monitor()
{
  delete timer_;
  timer_ = 0;
}

//----------------------------------------------------------------------

void Monitor::header () const
{
  print ("","==============================================");
  print ("","");
  print ("","  .oooooo.             oooo  oooo            ");
  print (""," d8P'  `Y8b            `888  `888            ");
  print ("","888           .ooooo.   888   888   .ooooo.  ");
  print ("","888          d88' `88b  888   888  d88' `88b ");
  print ("","888          888ooo888  888   888  888   888 ");
  print ("","`88b    ooo  888    .o  888   888  888   888 ");
  print (""," `Y8bood8P'  `Y8bod8P' o888o o888o `Y8bod8P' ");
  print ("","");
  print ("","A Parallel Adaptive Mesh Refinement Framework");
  print ("","");  
  print ("","  Laboratory for Computational Astrophysics");
  print ("","        San Diego Supercomputer Center");
  print ("","     University of California, San Diego");
  print ("","");  
  print ("","See 'LICENSE_CELLO' for software license information");
  print ("","");  

  // Get date text

  time_t rawtime;
  struct tm * t;
  time(&rawtime);
  t = localtime (&rawtime);
  const char * month[] = 
    {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

  print ("","BEGIN CELLO: %s %02d %02d:%02d:%02d",
	 month[t->tm_mon],
	 t->tm_mday,
	 t->tm_hour,
	 t->tm_min,
	 t->tm_sec);

  // Print all recognized configuration settings

  print ("Define","Simulation processors %d",CkNumPes());

  // Parallel type defines

  print ("Define","CELLO_ARCH %s",CELLO_ARCH);
  print ("Define","CELLO_PREC %s",CELLO_PREC);

  print ("Define","CC            %s",CELLO_CC);
  print ("Define","CFLAGS        %s",CELLO_CFLAGS);
  print ("Define","CPPDEFINES    %s",CELLO_CPPDEFINES);
  print ("Define","CPPPATH       %s",CELLO_CPPPATH);
  print ("Define","CXX           %s",CELLO_CXX);
  print ("Define","CXXFLAGS      %s",CELLO_CXXFLAGS);
  print ("Define","FORTRANFLAGS  %s",CELLO_FORTRANFLAGS);
  print ("Define","FORTRAN       %s",CELLO_FORTRAN);
  print ("Define","FORTRANLIBS   %s",CELLO_FORTRANLIBS);
  print ("Define","FORTRANPATH   %s",CELLO_FORTRANPATH);
  print ("Define","LIBPATH       %s",CELLO_LIBPATH);
  print ("Define","LINKFLAGS     %s",CELLO_LINKFLAGS);
  print ("Define","BUILD HOST    %s",CELLO_HOST);
  print ("Define","BUILD DIR     %s",CELLO_DIR);
  print ("Define","BUILD DATE    %s",CELLO_DATE);
  print ("Define","BUILD TIME    %s",CELLO_TIME);
  print ("Define","CHARM_VERSION %s",CELLO_CHARM_VERSION);

}

//----------------------------------------------------------------------

bool Monitor::is_active(const char * component) const throw ()
{
  if (! active_) return false;
  
  std::map<std::string,bool>::const_iterator it_active
    = group_active_.find(component);

  bool in_list = it_active != group_active_.end();

  return in_list ? it_active->second : group_default_;
}

//----------------------------------------------------------------------

void Monitor::write 
( FILE * fp, const char * component, const char * format,  ... ) const
{

  if (is_active(component)) {

    va_list fargs;

    // Process any input arguments

    char message[MONITOR_LENGTH+1];

    va_start(fargs,format);
    vsnprintf (message,MONITOR_LENGTH, format,fargs);
    va_end(fargs);

    write_ (fp, component,message);
  }
}

//----------------------------------------------------------------------

void Monitor::verbose 
( FILE * fp, const char * component, const char * format,  ... ) const
{

  if (verbose_ && is_active(component)) {

    va_list fargs;

    // Process any input arguments

    char message[MONITOR_LENGTH+1];
    va_start(fargs,format);
    vsnprintf (message,MONITOR_LENGTH, format,fargs);
    va_end(fargs);

    write_ (fp, component,message);
  }
}

//----------------------------------------------------------------------

void Monitor::write_ (FILE * fp, const char * component, const char * message) const
{

  // Get parallel process text

  char process[MONITOR_LENGTH] = "";
    
  sprintf (process,"%0d",CkMyPe());

  // Get time

  char time[10];

  snprintf (time,10,"%08.2f",timer_->value());

  // Print 

  if (fp == stdout) {
    PARALLEL_PRINTF 
      ("%s %s %s %s\n",
       process, time, component, message);
    fflush(stdout);
  } else {
    fprintf 
      (fp,"%s %s %s %s\n",
       process, time, component, message);
  }
}

//----------------------------------------------------------------------

void Monitor::write_verbatim
(
 FILE * fp,
 const char * component,
 const char * message
 ) const
{
  if (is_active(component)) {

    // Get parallel process text

    char buffer_process[MONITOR_LENGTH] = "";
    
    sprintf (buffer_process,"%0d",CkMyPe());

    // Get time

    char buffer_time[10];

    snprintf (buffer_time,10,"%08.2f",timer_->value());

    // Print 

    if (fp == stdout) {
      PARALLEL_PRINTF 
	("%s %s %s %s\n",
	 buffer_process, buffer_time, component, message);
      fflush(stdout);
    } else {
      fprintf 
	(fp,"%s %s %s %s\n",
	 buffer_process, buffer_time, component, message);
    }
  }

}

//----------------------------------------------------------------------

void Monitor::print (const char * component, const char * message, ...) const
{
  if (is_active(component)) {

    va_list fargs;

    // Process any input arguments

    char buffer_message[MONITOR_LENGTH+1];

    va_start(fargs,message);
    vsnprintf (buffer_message,MONITOR_LENGTH, message,fargs);
    va_end(fargs);

    write (stdout, component, buffer_message);
  }
}

//----------------------------------------------------------------------

void Monitor::print_verbatim (const char * component, const char * message) const
{
  write_verbatim (stdout, component, message);
}
