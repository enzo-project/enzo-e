// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
/// @bug      Probable memory leaks
/// @brief    Read in a parameter file and access parameter values

#include <stdlib.h>
#include <string.h>

#include "error.hpp"
#include "parameters.hpp"

Parameters::Parameters() 
  throw()
  :
  current_group_(""),
  current_subgroup_(""),
  parameter_map_(),
  parameter_list_()
  ///
{
}

Parameters::~Parameters()
///
{
  // Iterate over all parameters, deleting their values

  std::map<std::string,Param *>::iterator it_param;
  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {

    delete it_param->second;
  }
}

void Parameters::read ( FILE * file_pointer )
/// @param    file_pointer An opened input parameter file or stdin
{
  parameter_list_ = cello_parameters_read(file_pointer);

  struct param_type * node = parameter_list_ -> next; // skip sentinel
  struct param_type * prev = node;

  while (node->type != enum_parameter_sentinel) {

    std::string parameter_name = 
      std::string(node->group) + ":" +
      std::string(node->subgroup) + ":" +
      std::string(node->parameter);

    Param * param;

    param = new Param;

    param->set(node);

    parameter_map_     [parameter_name] = param;
    parameter_accessed_[parameter_name] = false;

    node = node->next;
    
    free (prev->group);
    free (prev->subgroup);
    free (prev->parameter);

    free (prev); // free not delete since allocated in parse.y

    prev = node;
    
  }
}

void Parameters::write ( FILE * file_pointer )
/// @param    file_pointer An opened output parameter file or stdout
{
  std::map<std::string,Param *>::iterator it_param;

  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {

    it_param->second->write(file_pointer, it_param->first);

  }
}

int Parameters::value_integer 
( std::string parameter,
  int         deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_integer()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_integer_ : deflt;
}


double Parameters::value_scalar 
( std::string parameter,
  double      deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return scalar (double) parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_scalar()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_scalar_ : deflt;
}

bool Parameters::value_logical 
( std::string parameter,
  bool        deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return logical parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_logical()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_logical_ : deflt;
}

std::string Parameters::value_string 
( std::string parameter,
  std::string deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_string()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_string_ : deflt;
}

void Parameters::evaluate_scalar 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated scalar parameters values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_scalar_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_scalar(param->value_expr_,n,result,x,y,z,t);
  } else {
    for (int i=0; i<n; i++) result[i] = deflt[i];
  }
}

void Parameters::evaluate_logical 
  (
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated logical parameters values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_logical_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_logical(param->value_expr_,n,result,x,y,z,t);
  } else {
    for (int i=0; i<n; i++) result[i] = deflt[i];
  }
}

int Parameters::list_length(std::string parameter)
/// @param   parameter Parameter name
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_list()) throw ExceptionParametersBadType();
  return (param != NULL) ? (param->value_list_)->size() : 0;
}

int Parameters::list_value_integer 
( int index,
  std::string parameter,
  int         deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the integer list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  if (param && ! param->is_integer()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_integer_ : deflt;
}

double Parameters::list_value_scalar 
( int index,
  std::string parameter,
  double      deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the scalar (double) list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return scalar (double) list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  if (param && ! param->is_scalar()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_scalar_ : deflt;
}

bool Parameters::list_value_logical 
( int index,
  std::string parameter,
  bool        deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the logical list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return logical list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  if (param && ! param->is_logical()) throw ExceptionParametersBadType();
  return (param != NULL) ? param->value_logical_ : deflt;
}

std::string Parameters::list_value_string 
( int index,
  std::string parameter,
  std::string deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the string list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string list parameter element value if it exists, deflt if not
{
  Param * param = list_element_ (parameter,index);
  if (param && ! param->is_string()) throw (ExceptionParametersBadType());
  return (param != NULL) ? param->value_string_ : deflt;
}

void Parameters::list_evaluate_scalar 
(
 int index,
 std::string parameter,
 int         n, 
 double    * result, 
 double    * deflt,
 double    * x, 
 double    * y, 
 double    * z, 
 double    * t
 ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated scalar expression list parameter element values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{

  Param * param = list_element_(parameter,index);
  if (param && ! param->is_scalar_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_scalar(param->value_expr_,n,result,x,y,z,t);
  } else {
    for (int i=0; i<n; i++) result[i] = deflt[i];
  }
}

void Parameters::list_evaluate_logical 
  (
   int index,
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated logical expression list parameter element values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{
  Param * param = list_element_(parameter,index);
  if (param && ! param->is_logical_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_logical(param->value_expr_,n,result,x,y,z,t);
  } else {
    for (int i=0; i<n; i++) result[i] = deflt[i];
  }
}


int Parameters::readline_ 
( FILE*  file_pointer, 
  char * buffer, 
  int    buffer_length ) 
  throw()
/// @param   file_pointer       the opened input file
/// @param   buffer             a character string buffer
/// @param   buffer_length      size of the buffer
/// @return  Whether the end of the file has been reached
{

  int i;

  // Clear the line buffer

  for (i=0; i<buffer_length; i++) {
    buffer[i]=0;
  }

  // Read the next line into the buffer

  int c = 0;
  for (i=0; c != EOF && c != '\n' && i < buffer_length-1; i++) {
    buffer[i] = c = fgetc(file_pointer);
  }

  // Back up i to last character read
  i--;

  // Check for buffer overrun

  if (i+1 >= buffer_length-1) {
    throw "Input file line too long";
  }

  // Convert the buffer into a C string

  if (buffer[i] == '\n') buffer[i] = '\0';

  // Return whether or not the end-of-file has been reached

  return (c != EOF);

}
