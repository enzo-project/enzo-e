// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoUserDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoUserDescr user-dependent class member functions

#include "cello.hpp"
#include "global.hpp"
#include "parameters.hpp"

#include "user.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoUserDescr::EnzoUserDescr(Global * global) throw ()
  : UserDescr(global),
    enzo_(new EnzoDescr(global))
{
  enzo_->read_parameters(global->parameters());
}
//----------------------------------------------------------------------

UserControl * 
EnzoUserDescr::create_user_control_ (std::string control_name) throw ()
{
  return new MethodEnzoControl(global_,enzo_);
}

//----------------------------------------------------------------------

UserTimestep * 
EnzoUserDescr::create_user_timestep_ ( std::string timestep_name ) throw ()
{
  return new MethodEnzoTimestep(enzo_);
}

//----------------------------------------------------------------------

UserMethod * 
EnzoUserDescr::create_user_method_ ( std::string method_name ) throw ()
/// @param method_name   Name of the method to create
{
  printf ("%s:%d '%s'\n",__FILE__,__LINE__,method_name.c_str());
  if (method_name == "ppm") {
    return new MethodEnzoPpm (global_,enzo_);
  } else {
    char buffer[80];
    sprintf (buffer,"Unknown method '%s'",method_name.c_str());
    global_->error()->error_ (__FILE__,__LINE__,"EnzoUserDescr::create_method",
				       buffer);
    return 0;
  }
}
