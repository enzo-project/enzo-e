// $Id: method_method.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 13 11:12:25 PDT 2009
/// @brief    Implements the Method base class

#include <string>

#include "method.hpp"
#include "monitor.hpp"

//----------------------------------------------------------------------

Method::Method() throw ()
///
{
  INCOMPLETE_MESSAGE("Method",""); 
}

//----------------------------------------------------------------------

void Method::add_argument_
(
 argument_type argument_type,
 std::string   argument_name,
 access_type   access_type
 ) throw()
/// @param         argument_type  Type of argument, field or particle
/// @param         argument_name  Name of the argument, e.g. "Density"
/// @param         access_type    Access type of the argument, e.g. read, write
{

  // Monitor output
  Monitor * monitor = Monitor::instance();
  char buffer[100];
  sprintf (buffer,"Method %s: adding %s", method_name_.c_str(), argument_name.c_str());
  monitor->print (buffer);

  // Add method argument information
  argument_types_.push_back(argument_type);
  argument_names_.push_back(argument_name);
  access_types_.push_back  (access_type);

}
