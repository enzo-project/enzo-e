#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      parameters.hpp
 * @brief     Read in a parameter file and access parameter values
 * @author    James Bordner
 * @date      Thu Jul  9 15:44:21 PDT 2009
 * @bug       
 * @note      
 *
 * Class Parameters to read in a parameter file and access
 * parameter values
 *
 * $Id$
 *
 *********************************************************************
 */

#include <map>
#include <string>
#include <stdio.h>

#include "cello.h"
#include "error_exception.hpp"
#include "type_parameter.h"
// Maximum allowed width of a line in a parameter file

#define MAX_PARAMETER_FILE_WIDTH 255

// Types for parameters


class Parameters {

/** 
 *********************************************************************
 *
 * @class     Parameters
 * @brief     Read in a parameter file and access parameter values
 * @ingroup   Parameters
 *
 * Read in a parameter file and access parameter values
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// This function creates an empty parameters object
  Parameters() throw();

  /// Read in parameters from a file

  void read (FILE * file_pointer) throw(ExceptionBadPointer);
  void read_bison (FILE * file_pointer) throw(ExceptionBadPointer);

  /// Return the string-valued parameter

  std::string value_string ( std::string parameter, 
			     std::string deflt = "") throw();

  /// Return the scalar-valued parameter

  Scalar value_scalar (std::string parameter, 
		       Scalar deflt = 0.0) throw();

  /// Return the integer-valued parameter

  int value_integer (std::string parameter, 
		     int deflt = 0) throw();

  /// Return the logical-valued parameter

  bool value_logical (std::string parameter, 
		      bool deflt = false) throw(ExceptionParametersBadType());

  /// Access parameters specific to a group or (group, subgroup)

  void set_group  (std::string group, std::string subgroup = "") throw ()
  { 
    current_group_ = group;
    current_subgroup_ = subgroup;
  };

  std::string get_group () throw ()
  { return current_group_; };
  std::string get_subgroup () throw ()
  { return current_subgroup_; };

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Add a parameter / value pair
  void add_parameter_ ( std::string parameter,  std::string value )   throw();

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// A private attribute

  std::map <std::string,std::string> values_;

  std::string current_group_;
  std::string current_subgroup_;
  
};

// Parser functions

extern "C" { 
  void cello_parameters_read(FILE *);
  void cello_parameters_print();
}

#endif /* PARAMETERS_HPP */

