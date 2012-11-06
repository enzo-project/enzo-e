// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"

#include "monitor.hpp" 

//----------------------------------------------------------------------
Monitor Monitor::instance_; // singleton design pattern)
//----------------------------------------------------------------------

Monitor::Monitor()
  : timer_(new Timer),
    active_(true),
    ip_(0)
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
  print ("","                James Bordner");
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

  // Print all compile-time configuration setings

  //  print ("Define","CELLO_ARCH = %s",CELLO_ARCH);
  //  print ("Define","CELLO_PREC = %s",CELLO_PREC);
  //  print ("Define","CELLO_TYPE = %s",CELLO_TYPE);

  // Print all recognized configuration settings

  // Parallel type defines

#ifdef CONFIG_USE_MPI
  print ("Define","CONFIG_USE_MPI");
#endif
#ifdef CONFIG_USE_CHARM
  print ("Define","CONFIG_USE_CHARM");
#endif

  // Precision defines

#ifdef CONFIG_PRECISION_SINGLE
  print ("Define","CONFIG_PRECISION_SINGLE");
#endif
#ifdef CONFIG_PRECISION_DOUBLE
  print ("Define","CONFIG_PRECISION_DOUBLE");
#endif

  // Performance defines

#ifdef CONFIG_USE_MEMORY
  print ("Define","CONFIG_USE_MEMORY");
#endif
#ifdef CONFIG_USE_PROJECTIONS
  print ("Define","CONFIG_USE_PROJECTIONS");
#endif
#ifdef CONFIG_USE_PERFORMANCE
  print ("Define","CONFIG_USE_PERFORMANCE");
#endif
#ifdef CONFIG_USE_PAPI
  print ("Define","CONFIG_USE_PAPI");
#endif

  // Debugging defines

#ifdef CELLO_TRACE
  print ("Define","CELLO_TRACE");
#endif
#ifdef CELLO_DEBUG
  print ("Define","CELLO_DEBUG");
#endif
#ifdef CELLO_DEBUG_VERBOSE
  print ("Define","CELLO_DEBUG_VERBOSE");
#endif

  // Library defines

#ifdef H5_USE_16_API
  print ("Define","H5_USE_16_API");
#endif
#ifdef NO_FREETYPE
  print ("Define","NO_FREETYPE");
#endif

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
(
 FILE * fp,
 const char * component,
 const char * message,
  ...
 ) const
{

  if (is_active(component)) {

    va_list fargs;

    // Process any input arguments

    char buffer_message[MONITOR_LENGTH+1];

    va_start(fargs,message);
    vsnprintf (buffer_message,MONITOR_LENGTH, message,fargs);
    va_end(fargs);
    
    // Get parallel process text

    char buffer_process[MONITOR_LENGTH] = "";
    
    sprintf (buffer_process,"%0d",ip_);

    // Get time

    char buffer_time[10];

    snprintf (buffer_time,10,"%08.2f",timer_->value());

    // get Component if any
    char buffer_component[20];

    if (strlen(component)>0) {
      snprintf (buffer_component,20," %-11s ",component);
    } else {
      buffer_component[0] = 0;
    }

    // Print 

    if (fp == stdout) {
      PARALLEL_PRINTF 
	("%s %s %s %s\n",
	 buffer_process, buffer_time, buffer_component, buffer_message);
    } else {
      fprintf 
	(fp,"%s %s %s %s\n",
	 buffer_process, buffer_time, buffer_component, buffer_message);
    }

#ifdef CELLO_DEBUG
    // Write DEBUG messagess to file

    if (strcmp(component,"DEBUG")==0) {
      char file[20];
      sprintf (file,"out.debug.%d",ip_);
      FILE * fdebug = fopen (file,"a");
      fprintf (fdebug,"%s %s %s %s\n",
	 buffer_process, buffer_time, buffer_component, buffer_message);
      fclose(fdebug);
    }

#endif

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
    
    sprintf (buffer_process,"%0d",ip_);

    // Get time

    char buffer_time[10];

    snprintf (buffer_time,10,"%08.2f",timer_->value());

    // get Component if any
    char buffer_component[20];

    if (strlen(component)>0) {
      snprintf (buffer_component,20," %-11s ",component);
    } else {
      buffer_component[0] = 0;
    }

    // Print 

    if (fp == stdout) {
      PARALLEL_PRINTF 
	("%s %s %s %s\n",
	 buffer_process, buffer_time, buffer_component, message);
    } else {
      fprintf 
	(fp,"%s %s %s %s\n",
	 buffer_process, buffer_time, buffer_component, message);
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
