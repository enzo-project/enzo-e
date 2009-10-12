//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      parameters.cpp
 * @brief     Read in a parameter file and access parameter values
 * @author    James Bordner
 * @date      Thu Jul  9 15:38:43 PDT 2009
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Class Parameters to read in a parameter file and access
 *    parameter values
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    parameters.hpp
 *
 * PUBLIC FUNCTIONS
 *  
 *    Parameters::Parameters()
 *    Parameters::~Parameters()
 *    Parameters::read()
 *    Parameters::value_string()
 *    Parameters::value_scalar()
 *    Parameters::value_integer()
 *    Parameters::value_logical()
 *
 * PRIVATE FUCTIONS
 *  
 * $Id$
 *
 *********************************************************************
 */

#include <stdlib.h>
#include <string.h>

#include "error.hpp"
#include "parameters.hpp"

/**
 *********************************************************************
 *
 * @param   
 * @return  
 *
 * This function creates an empty Parameters object
 *
 *********************************************************************
 */

Parameters::Parameters() 
  throw()
  :
  current_group_("testing"),
  current_subgroup_("testing"),
  parameter_list_(NULL)
{
}

/**
 *********************************************************************
 *
 * @param   
 * @return  
 *
 * This function deletes a Parameters object
 *
 *********************************************************************
 */

Parameters::~Parameters()
{
  std::map<std::string,Param *>::iterator it_param;
  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {
    it_param->second->dealloc();
    delete it_param->second;
  }
}

/**
 *********************************************************************
 *
 * @param   file_pointer: An opened input parameter file or stdin
 * @return  There is no return value
 *
 * This function reads in parameter-value key pairs, one per line
 *
 *********************************************************************
 */


void
Parameters::read ( FILE * file_pointer ) throw(ExceptionBadPointer)

{
  parameter_list_ = cello_parameters_read(file_pointer);

  struct param_type * node = parameter_list_ -> next; // skip sentinel
  struct param_type * prev = node;

  while (node->type != enum_parameter_sentinel) {

    std::string parameter_name = 
      std::string(node->group) + ":" +
      std::string(node->subgroup) + ":" +
      std::string(node->parameter);

    // monitor: debug
    //    if (debug) printf ("parameter %s\n",parameter_name.c_str()); fflush(stdout);

    enum_parameter type = node->type;
    
    Param * param;

    switch (type) {
    case enum_parameter_integer:
      param = new Param;
      param->set_integer_(node->integer_value);
      break;
    case enum_parameter_scalar:
      param = new Param;
      param->set_scalar_(node->scalar_value);
      break;
    case enum_parameter_string:
      param = new Param;
      param->set_string_(node->string_value);
      break;
    case enum_parameter_logical:
      param = new Param;
      param->set_logical_(node->logical_value);
      break;
    case enum_parameter_list:
      param = new Param;
      param->set_list_(node->list_value);
      break;
    case enum_parameter_scalar_expr:
      param = new Param;
      param->set_scalar_expr_(node->op_value);
      break;
    case enum_parameter_logical_expr:
      param = new Param;
      param->set_logical_expr_(node->op_value);
      break;
    }

    parameter_map_[parameter_name] = param;

    node = node->next;
    
    free (prev->group);
    free (prev->subgroup);
    free (prev->parameter);

    free (prev); // free not delete since allocated in parse.y

    prev = node;
    
  }
}

/**
 *********************************************************************
 *
 * @param   parameter
 * @return  Return the string value of the parameter if it exists, 
 *          or deflt if not
 *
 * Return the string value of the parameter if it exists, or deflt if
 * not.
 *
 *********************************************************************
 */

std::string 
Parameters::value_string 
( std::string parameter,
  std::string deflt ) throw()
{
  Param * param = parameter_(parameter);
  return (param != NULL) ? param->value_string_ : deflt;
}

/**
 *********************************************************************
 *
 * @param   parameter
 * @return  Return the scalar value of the parameter if it exists, 
 *          or deflt if not
 *
 * Return the scalar value of the parameter if it exists, or deflt if
 * not.
 *
 *********************************************************************
 */

Scalar
Parameters::value_scalar 
( std::string parameter,
  Scalar      deflt ) throw()
{
  Param * param = parameter_(parameter);
  return (param != NULL) ? param->value_scalar_ : deflt;
}

/**
 *********************************************************************
 *
 * @param   parameter
 * @return  Return the integer value of the parameter if it exists, 
 *          or deflt if not
 *
 * Return the integer value of the parameter if it exists, or deflt if
 * not.
 *
 *********************************************************************
 */

int
Parameters::value_integer 
( std::string parameter,
  int         deflt ) throw()
{
  Param * param = parameter_(parameter);
  return (param != NULL) ? param->value_integer_ : deflt;
}

/**
 *********************************************************************
 *
 * @param   parameter
 * @return  Return the logical value of the parameter if it exists, 
 *          or deflt if not
 *
 * Return the logical value of the parameter if it exists, or deflt if
 * not.
 *
 *********************************************************************
 */

bool
Parameters::value_logical 
( std::string parameter,
  bool        deflt ) throw(ExceptionParametersBadType())
{
  Param * param = parameter_(parameter);
  return (param != NULL) ? param->value_logical_ : deflt;
}

/**
 *********************************************************************
 *
 * @param   parameter
 *
 * @return  Return the scalar value of the expression parameter if it
 *          exists, or deflt if not
 *
 * Return the scalar value of the expression parameter if it exists,
 * or deflt if not.
 *
 *********************************************************************
 */

void Parameters::evaluate_scalar 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
{
  Param * param = parameter_(parameter);
  if (param != NULL) {
    param->evaluate_scalar(n,result,x,y,z,t);
  } else {
    for (int i=0; i<n; i++) result[i] = deflt[i];
  }
}

//======================================================================
// PRIVATE FUNCTIONS
//======================================================================

/**
 *********************************************************************
 *
 * @param   file_pointer       the opened input file
 * @param   buffer             a character string buffer
 * @param   buffer_length      size of the buffer
 *
 * @return  Whether the end of the file has been reached
 *
 * Read in the next line of the input file into the buffer
 *
 *********************************************************************
 */

int 
Parameters::readline_ 
( FILE*  file_pointer, 
  char * buffer, 
  int    buffer_length ) 
  throw()

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

void    
Parameters::add_parameter_ 
( std::string parameter,
  std::string value ) 
  throw()

{
//     if (parameter=="" || parameter=="#" || parameter=="//") return;
//   values_.insert( std::pair<std::string,std::string>(parameter,value) );
}
