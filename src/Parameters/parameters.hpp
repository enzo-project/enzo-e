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
#include <stdlib.h>

#include "cello.h"
#include "error.hpp"

#include "parse.h"
#include "param.hpp"

// Maximum allowed width of a line in a parameter file

#define MAX_PARAMETER_FILE_WIDTH 255

// Functions for deleting parse.* structs

// Types for parameters

class Parameters {


  //   friend class Parameter;

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

  /// Create an empty Parameters object

  Parameters() throw();

  /// Deletes a Parameters object

  ~Parameters();

  /// Read in parameters from a file

  void read (FILE * file_pointer) throw(ExceptionBadPointer);

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

  int value_list_length (std::string parameter) {
  }

  void evaluate_scalar 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t);

  /// Access parameters specific to a group or (group, subgroup)

  void set_group  (std::string group, std::string subgroup = "") throw ()
  { 
    current_group_ = group;
  };

  void set_subgroup  (std::string subgroup) throw ()
  {
    current_subgroup_ = subgroup;
  }

  std::string get_group () throw ()
  { return current_group_; };
  std::string get_subgroup () throw ()
  { return current_subgroup_; };

private:

  //-------------------------------------------------------------------
  // PRIVATE FUNCTIONS
  //-------------------------------------------------------------------

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Add a parameter / value pair
  void add_parameter_ ( std::string parameter,  std::string value )   throw();

  /// Return the Param pointer for the specified parameter
  Param * parameter_ (std::string parameter)
  {
    std::string p = current_group_ + ":" + current_subgroup_ + ":" + parameter;
    return parameter_map_[p];
  };


private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// A private attribute

  std::string current_group_;
  std::string current_subgroup_;

  std::map<std::string, Param *> parameter_map_;
  struct param_type * parameter_list_;
};

// Parser functions

extern "C" { 
  struct param_type * cello_parameters_read(FILE *);
  void cello_parameters_print();
}

#endif /* PARAMETERS_HPP */

