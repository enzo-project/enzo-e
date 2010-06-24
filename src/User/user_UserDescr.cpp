// $Id: user_UserDescr.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_UserDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of UserDescr user-dependent class member functions

#include "user.hpp"
#include "error.hpp"

//----------------------------------------------------------------------

UserControl * UserDescr::create_user_control_ (std::string control_name)
{
  return new MethodEnzoControl;
}

//----------------------------------------------------------------------

UserTimestep * UserDescr::create_user_timestep_ (std::string timestep_name)
{
  return new MethodEnzoTimestep;
}

//----------------------------------------------------------------------

UserMethod * UserDescr::create_user_method_ (std::string method_name)
/// @param method_name   Name of the method to create
{
  if (method_name == "ppm") {
    return new MethodEnzoPpm;
  } else {
    char buffer[80];
    sprintf (buffer,"Unknown method %s",method_name.c_str());
    WARNING_MESSAGE ("UserDescr::create_method",
		     buffer);
    return 0;
  }
}
