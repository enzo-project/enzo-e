// $Id: user_MethodDescr.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_MethodDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of MethodDescr user-dependent class member functions

#include "method.hpp"
#include "user.hpp"
#include "error.hpp"

//----------------------------------------------------------------------

MethodControl * MethodDescr::create_method_control_ (std::string method_name)
{
  return new MethodEnzoControl;
}

//----------------------------------------------------------------------

MethodTimestep * MethodDescr::create_method_timestep_ (std::string method_name)
{
  return new MethodEnzoTimestep;
}

//----------------------------------------------------------------------

MethodUser * MethodDescr::create_method_user_ (std::string method_name)
/// @param method_name   Name of the method to create
{
  if (method_name == "ppm") {
    return new MethodEnzoPpm;
  } else {
    char buffer[80];
    sprintf (buffer,"Unknown method %s",method_name.c_str());
    WARNING_MESSAGE ("MethodDescr::create_method",
		     buffer);
    return 0;
  }
}
