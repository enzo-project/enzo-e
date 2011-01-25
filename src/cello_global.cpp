// $Id: cello_global.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     cello_global.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file cello_global.cpp
///
/// Constant static global variables 


#include "cello.hpp"

#include "error.hpp"
#include "monitor.hpp"

//----------------------------------------------------------------------

const char * component_name [] = {
  // !!! EDIT component_type AND component_name TOGETHER !!!
  "undefined",
  "enzop",
  "disk",
  "distribute",
  "error",
  "field",
  "memory",
  "mesh",
  "method",
  "monitor",
  "parallel",
  "parameters",
  "particles",
  "performance",
  "portal",
  "schedule",
  "simulation",
  "task"
};

// //----------------------------------------------------------------------

// Global::Global() throw()
// {
//   error_      = new Error;
//   monitor_    = new Monitor;
// //   parameters_ = new Parameters(monitor_);
//   //    memory_     = new Memory;
// }

// //----------------------------------------------------------------------

// Global::~Global() throw()
// {
//   delete error_;
// //   delete parameters_;
//   delete monitor_;
//   //    delete [] memory_;
// }
