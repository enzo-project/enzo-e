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

  /// Delete a Parameters object

  ~Parameters();

  /// Read in parameters from a file

  void read (FILE * file_pointer) 
    throw(ExceptionBadPointer);

  /// Return the integer-valued parameter

  int value_integer (std::string , int deflt = 0) 
    throw(ExceptionParametersBadType);

  /// Return the scalar-valued parameter

  double value_scalar (std::string, double deflt = 0.0) 
    throw(ExceptionParametersBadType);

  /// Return the logical-valued parameter

  bool value_logical (std::string , bool deflt = false) 
    throw(ExceptionParametersBadType);

  /// Return the string-valued parameter

  std::string value_string ( std::string , std::string deflt = "") 
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued parameter expression

  void evaluate_scalar 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
    throw(ExceptionParametersBadType);

  /// Evaluate the logical-valued parameter expression

  void evaluate_logical 
  (
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
    throw(ExceptionParametersBadType);

  /// Return the length of the list parameter

  int list_length (std::string parameter);

  /// Access an integer list element
  int list_value_integer (int , std::string , int deflt = 0)    
    throw(ExceptionParametersBadType);

  /// Access a scalar list element
  double list_value_scalar (int , std::string , double deflt = 0.0)    
    throw(ExceptionParametersBadType);

  /// Access a logical list element
  bool list_value_logical (int ,std::string , bool deflt = false)    
    throw(ExceptionParametersBadType);

  /// Access a string list element
  std::string list_value_string (int ,std::string , std::string deflt = "")    
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued list element expression
  void list_evaluate_scalar 
  (
   int ,
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t
   )    
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued list element expression
  void list_evaluate_logical
  (
   int ,
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t
   )    
    throw(ExceptionParametersBadType);

  /// Set the current group

  void set_group  (std::string group, std::string subgroup = "") throw ()
  { 
    current_group_    = group;
    current_subgroup_ = "";
  };

  /// Set the current subgroup

  void set_subgroup  (std::string subgroup) throw ()
  {
    current_subgroup_ = subgroup;
  }

  /// Get the current group

  std::string get_group () throw ()
  { return current_group_; };

  /// Get the current subgroup

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

  /// Return the Param pointer for the specified list parameter element

  Param * list_element_ (std::string parameter, int index)
  {
    Param * list = parameter_(parameter);
    Param * param = NULL;
    if (list != NULL) {
      param = (*(list->value_list_))[index];
    }
    return param;
  }


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

