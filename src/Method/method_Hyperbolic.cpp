// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Hyperbolic.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 13 11:12:25 PDT 2009
/// @brief    Implements the Hyperbolic base class

#include "cello.hpp"

#include "method.hpp"

//----------------------------------------------------------------------

void Hyperbolic::add_argument_
(
 argument_enum argument,
 std::string   argument_name,
 access_enum   access_type,
 DataDescr   * data_descr
 ) throw()
/// @param    argument       Type of argument, field or particle
/// @param    argument_name  Name of the argument, e.g. "Density"
/// @param    access_type    Access type of the argument, e.g. read, write
/// @param    data_descr     The data descriptor
{

  // Monitor output

  char buffer[100];
  sprintf (buffer,"Hyperbolic %s: adding %s", 
	   name().c_str(), 
	   argument_name.c_str());
  monitor_->print (buffer);

  // Add method argument information
  argument_types_.push_back(argument);
  argument_names_.push_back(argument_name);
  access_types_.push_back  (access_type);

  // If data_descr is passed in (default = 0), then verify that the argument
  // is defined
  if (data_descr) {
    char buffer [ ERROR_LENGTH ];
    switch (argument) {
    case argument_field:
      sprintf (buffer, 
	       "Required Field %s is not defined in the field descriptor",
	       argument_name.c_str());
      ASSERT("Hyperbolic::add_argument",
	     buffer, data_descr->field_descr()->is_field(argument_name));
      break;
    case argument_particle:
    case argument_unknown:
      break;
    }
  }

}
