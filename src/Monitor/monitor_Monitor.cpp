// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @todo      simplify image call
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"

#include "monitor.hpp" 

// Monitor * Monitor::instance_ = 0; // (singleton design pattern)

//----------------------------------------------------------------------

void Monitor::header ()
{
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
  print ("James Bordner");
  print ("Laboratory for Computational Astrophysics");
  print ("San Diego Supercomputer Center");
  print ("University of California, San Diego");
  print ("");  
}

//----------------------------------------------------------------------

void Monitor::config () throw()
{
  print ("Configuration");
  print ("-------------");
#ifdef CONFIG_PRECISION_DOUBLE
  print ("CONFIG_PRECISION_DOUBLE");
#endif
#ifdef CONFIG_PRECISION_QUADRUPLE
  print ("CONFIG_PRECISION_QUADRUPLE");
#endif
#ifdef CONFIG_PRECISION_SINGLE
  print ("CONFIG_PRECISION_SINGLE");
#endif
#ifdef CONFIG_USE_CHARM
  print ("CONFIG_USE_CHARM");
#endif
#ifdef CONFIG_USE_HDF5
  print ("CONFIG_USE_HDF5");
#endif
#ifdef CONFIG_USE_MEMORY
  print ("CONFIG_USE_MEMORY");
#endif
#ifdef CONFIG_USE_MPI
  print ("CONFIG_USE_MPI");
#endif
#ifdef CONFIG_USE_PAPI
  print ("CONFIG_USE_PAPI");
#endif
  print ("-------------");
}

