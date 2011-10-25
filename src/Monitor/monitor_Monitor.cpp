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
    active_(true)
{ 
  timer_->start();
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

  char c_single = ' ';
  char c_double = ' ';
  char c_quad   = ' ';
  char c_charm  = ' ';
  char c_mpi    = ' ';
  char c_papi   = ' ';

#ifdef CONFIG_PRECISION_SINGLE
  c_single = '*';
#endif
#ifdef CONFIG_PRECISION_DOUBLE
  c_double = '*';
#endif
#ifdef CONFIG_PRECISION_QUAD
  c_quad = '*';
#endif
#ifdef CONFIG_USE_CHARM
  c_charm = '*';
#endif
#ifdef CONFIG_USE_MPI
  c_mpi = '*';
#endif
#ifdef CONFIG_USE_PAPI
  c_papi = '*';
#endif

  char s_single[80];
  char s_double[80];
  char s_quad  [80];
  char s_charm [80];
  char s_mpi   [80];
  char s_papi  [80];

  sprintf (s_single,"(%c) CONFIG_PRECISION_SINGLE",c_single);
  sprintf (s_double,"(%c) CONFIG_PRECISION_DOUBLE",c_double);
  sprintf (s_quad,  "(%c) CONFIG_PRECISION_QUAD",  c_quad);
  sprintf (s_charm, "(%c) CONFIG_USE_CHARM",       c_charm);
  sprintf (s_mpi,   "(%c) CONFIG_USE_MPI",         c_mpi);
  sprintf (s_papi,  "[%c] CONFIG_USE_PAPI",        c_papi);

  print ("","==============================================");
  print ("",s_single);
  print ("",s_double);
  print ("",s_quad);
  print ("","");
  print ("",s_charm);
  print ("",s_mpi);
  print ("","");
  print ("",s_papi);
  print ("","==============================================");

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
  if (active_) {

    va_list fargs;

    // Process any input arguments

    char buffer_message[MONITOR_LENGTH+1];

    va_start(fargs,message);
    vsnprintf (buffer_message,MONITOR_LENGTH, message,fargs);
    va_end(fargs);
    
    // Get parallel process text

    char buffer_process[MONITOR_LENGTH] = "";

#if defined(CONFIG_USE_CHARM)
    sprintf (buffer_process,"%0d",CkMyPe());
#elif defined(CONFIG_USE_MPI)
    sprintf (buffer_process,"%0d",Mpi::rank());
#endif

    // Get time

    char buffer_time[10];

    snprintf (buffer_time,10,"%08.2f",timer_->value());

    // get Component if any
    char buffer_component[20];

    if (strlen(component)>0) {
      snprintf (buffer_component,20," %-10s ",component);
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
  }

}

//----------------------------------------------------------------------

void Monitor::print (const char * component, const char * message, ...) const
{
  if (active_) {

    va_list fargs;

    // Process any input arguments

    char buffer_message[MONITOR_LENGTH+1];

    va_start(fargs,message);
    vsnprintf (buffer_message,MONITOR_LENGTH, message,fargs);
    va_end(fargs);

    write (stdout, component, buffer_message);
  }
}
