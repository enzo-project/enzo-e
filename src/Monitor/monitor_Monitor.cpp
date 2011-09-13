// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @todo      simplify image call
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"

#include "monitor.hpp" 

//----------------------------------------------------------------------
Monitor * Monitor::instance_ = 0; // (singleton design pattern)
//----------------------------------------------------------------------

Monitor::Monitor()
  : timer_(new Timer),
    active_(true),
    image_(0),
    image_size_x_(0),
    image_size_y_(0),
    png_(0)
{ 
  timer_->start();

  map_r_.resize(2);
  map_g_.resize(2);
  map_b_.resize(2);
  map_r_[0] = 0.0;
  map_g_[0] = 0.0;
  map_b_[0] = 0.0;
  map_r_[1] = 1.0;
  map_g_[1] = 1.0;
  map_b_[1] = 1.0;
}

//----------------------------------------------------------------------

Monitor::~Monitor()
{
  delete timer_;
  timer_ = 0;
  delete instance_;
  instance_ = 0;
}

//----------------------------------------------------------------------

void Monitor::header () const
{
  print ("==============================================");
  print ("");
  print ("  .oooooo.             oooo  oooo            ");
  print (" d8P'  `Y8b            `888  `888            ");
  print ("888           .ooooo.   888   888   .ooooo.  ");
  print ("888          d88' `88b  888   888  d88' `88b ");
  print ("888          888ooo888  888   888  888   888 ");
  print ("`88b    ooo  888    .o  888   888  888   888 ");
  print (" `Y8bood8P'  `Y8bod8P' o888o o888o `Y8bod8P' ");
  print ("");
  print ("A Parallel Adaptive Mesh Refinement Framework");
  print ("");  
  print ("                James Bordner");
  print ("  Laboratory for Computational Astrophysics");
  print ("        San Diego Supercomputer Center");
  print ("     University of California, San Diego");
  print ("");  

  // Get date text

  char buffer_date[MONITOR_LENGTH];

  time_t rawtime;
  struct tm * t;
  time(&rawtime);
  t = localtime (&rawtime);
  const char * month[] = 
    {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

  print ("BEGIN CELLO: %s %02d %02d:%02d:%02d",
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

  print ("==============================================");
  print (s_single);
  print (s_double);
  print (s_quad);
  print ("");
  print (s_charm);
  print (s_mpi);
  print ("");
  print (s_papi);
  print ("==============================================");

}

//----------------------------------------------------------------------

void Monitor::print (const char * message, ...) const
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

    // Print 

    PARALLEL_PRINTF ("%s %s %s\n",
		     buffer_process,
		     buffer_time,
		     buffer_message);
  }
}

//----------------------------------------------------------------------

void Monitor::image_set_map 
(int n, double * map_r, double * map_g, double * map_b) throw()
{
  map_r_.resize(n);
  map_g_.resize(n);
  map_b_.resize(n);

  for (int i=0; i<n; i++) {
    map_r_[i] = map_r[i];
    map_g_[i] = map_g[i];
    map_b_[i] = map_b[i];
  }
}

//----------------------------------------------------------------------

void Monitor::image_open (std::string filename, int mx, int my)
{
  png_ = new pngwriter(mx,my,0,filename.c_str());

  image_size_x_ = mx;
  image_size_y_ = my;

  image_ = new double [mx*my];

  for (int i=0; i<mx*my; i++) image_[i] = 0.0;
}

//----------------------------------------------------------------------

void Monitor::image_close (double min, double max)
{

  // simplified variable names

  int mx = image_size_x_;
  int my = image_size_y_;
  int m  = mx*my;

  // Adjust min and max bounds if needed

  for (int i=0; i<m; i++) {
    min = MIN(min,image_[i]);
    max = MAX(max,image_[i]);
  }

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_[i];

      // map v to lower colormap index
      size_t k = (map_r_.size() - 1)*(value - min) / (max-min);

      // prevent k == map_.size()-1, which happens if value == max
      if (k > map_r_.size() - 2) k = map_r_.size()-2;

      // linear interpolate colormap values
      double lo = min +  k   *(max-min)/(map_r_.size()-1);
      double hi = min + (k+1)*(max-min)/(map_r_.size()-1);

      // should be in bounds, but force if not due to rounding error
      if (value < lo) value = lo;
      if (value > hi) value = hi;

      // interpolate colormap

      double ratio = (value - lo) / (hi-lo);

      double r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
      double g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
      double b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];

      png_->plot(ix+1, iy+1, r,g,b);
    }
  }      

  png_->close();

  delete [] image_;
  image_ = 0;
}
