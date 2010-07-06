// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @todo      simplify image call
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"
#include "parallel.hpp"
#include "monitor.hpp" 
#include "error.hpp" 

// Monitor * Monitor::instance_ = 0; // (singleton design pattern)

void Monitor::header ()
{
  //    print ("");
  //    print ("     The Laboratory for Computational Astrophysics proudly presents:");
  print ("");
  print ("    =================================================================");
  print ("");
  print ("    oooooooooooo                                          ooooooooo.   ");
  print ("    `888'     `8                                          `888   `Y88. ");
  print ("     888         ooo. .oo.     oooooooo  .ooooo.           888   .d88' ");
  print ("     888oooo8    `888P\"Y88b   d'\"\"7d8P  d88' `88b          888ooo88P'  ");
  print ("     888    \"     888   888     .d8P'   888   888 8888888  888         ");
  print ("     888       o  888   888   .d8P'  .P 888   888          888         ");
  print ("    o888ooooood8 o888o o888o d8888888P  `Y8bod8P'         o888o        ");
  print ("");
  print ("    =================================================================");
  print ("              E N Z O : T H E   N E X T  G E N E R A T I O N");
  print ("    =================================================================");
  print ("");
}

