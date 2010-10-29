// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
/// @bug      Probable memory leaks
/// @brief    Read in a parameter file and access parameter values

#include <stdlib.h>
#include <string.h>

#include "cello.hpp"

#include "error.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

// Parameters Parameters::instance_; // (singleton design pattern)

//----------------------------------------------------------------------

Parameters::Parameters(Monitor * monitor) 
  throw()
  :  monitor_(monitor),
     current_group_(""),
     current_subgroup_(""),
     parameter_map_(),
     parameter_tree_(new ParamNode("Parameters"))
     ///
{
}

//----------------------------------------------------------------------

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
  delete parameter_tree_;
}

//----------------------------------------------------------------------

void Parameters::read ( FILE * file_pointer )
/// @param    file_pointer An opened input parameter file or stdin
{
  if ( file_pointer == NULL) {
    ERROR_MESSAGE("Parameters::read","file pointer is NULL");
  }
  
  struct param_struct * parameter_list = cello_parameters_read(file_pointer);

  struct param_struct * node = parameter_list -> next; // skip sentinel
  struct param_struct * prev = node;

  while (node->type != enum_parameter_sentinel) {

    Param * param = new Param;

    param->set(node);

    new_param_(node->group,node->subgroup,node->parameter,param);

    node = node->next;
    
    free (prev->group);
    free (prev->subgroup);
    free (prev->parameter);

    free (prev); // free not delete since allocated in parse.y

    prev = node;
    
  }
}

//----------------------------------------------------------------------

void Parameters::write ( FILE * file_pointer )
/// @param    file_pointer An opened output parameter file or stdout
{
  std::map<std::string,Param *>::iterator it_param;

  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {

    if (it_param->second) {
      it_param->second->write(file_pointer, it_param->first);
    } else {
      char message [ ERROR_MESSAGE_LENGTH ];
      sprintf (message, 
	       "uninitialized parameter %s accessed\n",
	       it_param->first.c_str());
      WARNING_MESSAGE("Parameters::write",message);
    }

  }
}

//----------------------------------------------------------------------

void Parameters::value 
(
 std::string    parameter, 
 parameter_type type, 
 void *         value, 
 void *         deflt
)
/// @param   parameter Parameter name
/// @param   type      Parameter type
/// @param   value     parameter value if it exists, deflt if not
/// @param   deflt     Default parameter value
{
  Param * param = parameter_(parameter);
  // Check type of parameter if it's set
  if (param) {
    bool correct_type;
    switch (type) {
    case parameter_integer:      correct_type = param->is_integer();      break;
    case parameter_scalar:       correct_type = param->is_scalar();       break;
    case parameter_logical:      correct_type = param->is_logical();      break;
    case parameter_string:       correct_type = param->is_string();       break;
    case parameter_list:
    case parameter_scalar_expr:
    case parameter_logical_expr:
      throw ExceptionParametersBadType();
      break;
    }
  }
  if (param != NULL) {
    switch (type) {
    case parameter_integer:
      *((int *)   value) = param->get_integer();
      break;
    case parameter_scalar:
      *((double *)value) = param->get_scalar();
      break;
    case parameter_logical:
      *((bool *)  value) = param->get_logical();
      break;
    case parameter_string:
      *((const char **) value) = param->get_string();
      break;
    case parameter_list:         
    case parameter_scalar_expr:  
    case parameter_logical_expr: 
      ERROR_MESSAGE ("Parameters::value","Switch case should never be reached");
      break;
    }
  } else if (deflt) {
    switch (type) {
    case parameter_integer:
      *((int *)value) = *((int *)deflt);
      break;
    case parameter_scalar:
      *((double *)value) = *((double *)deflt);
      break;
    case parameter_logical:
      *((bool *)value) = *((bool *)deflt);
      break;
    case parameter_string:
      *((const char **)value) = *((const char **)deflt);
      break;
    case parameter_list:         
    case parameter_scalar_expr:  
    case parameter_logical_expr: 
      ERROR_MESSAGE ("Parameters::value","Switch case should never be reached");
      break;
    }
  } else {
    char buffer[ERROR_MESSAGE_LENGTH];
    sprintf (buffer,"Required parameter %s not defined",parameter.c_str());
    ERROR_MESSAGE ("Parameters::value",buffer);
  }
  monitor_log(parameter);
}
//----------------------------------------------------------------------

int Parameters::value_integer 
( std::string parameter,
  int         deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_integer()) throw ExceptionParametersBadType();
  monitor_log(parameter);
  return (param != NULL) ? param->get_integer() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_integer 
( std::string parameter,
  int         value ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);
  if (param) {
    if (! param->is_integer()) throw ExceptionParametersBadType();
  } else {
    param = new Param;
    new_param_ (current_group_,current_subgroup_,parameter,param);
  }
  param->set_integer_(value);
  monitor_log(parameter);
}

//----------------------------------------------------------------------

double Parameters::value_scalar 
( std::string parameter,
  double      deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return scalar (double) parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_scalar()) throw ExceptionParametersBadType();
  monitor_log(parameter);
  return (param != NULL) ? param->get_scalar() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_scalar
( std::string parameter,
  double      value ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);
  if (param) {
    if (! param->is_scalar()) throw ExceptionParametersBadType();
  } else {
    param = new Param;
    new_param_ (current_group_,current_subgroup_,parameter,param);
  }
  param->set_scalar_(value);
  monitor_log(parameter);
}

//----------------------------------------------------------------------

bool Parameters::value_logical 
( std::string parameter,
  bool        deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return logical parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_logical()) throw ExceptionParametersBadType();
  monitor_log(parameter);
  return (param != NULL) ? param->get_logical() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_logical
( std::string parameter,
  bool        value ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);
  if (param) {
    if (! param->is_logical()) throw ExceptionParametersBadType();
  } else {
    param = new Param;
    new_param_ (current_group_,current_subgroup_,parameter,param);
  }
  param->set_logical_(value);
  monitor_log(parameter);
}

//----------------------------------------------------------------------

const char * Parameters::value_string 
( std::string  parameter,
  const char * deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_string()) throw ExceptionParametersBadType();
  monitor_log(parameter);
  return (param != NULL) ? param->get_string() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_string
( std::string  parameter,
  const char * value ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);
  if (param) {
    if (! param->is_string()) throw ExceptionParametersBadType();
  } else {
    param = new Param;
    new_param_ (current_group_,current_subgroup_,parameter,param);
  }
  param->set_string_(strdup(value));
  monitor_log(parameter);
}

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

int Parameters::list_length(std::string parameter)
/// @param   parameter Parameter name
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_list()) throw ExceptionParametersBadType();
  return (param != NULL) ? (param->value_list_)->size() : 0;
}

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

const char * Parameters::list_value_string 
( int index,
  std::string   parameter,
  const char *  deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the string list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string list parameter element value if it exists, deflt if not
{
  Param * param = list_element_ (parameter,index);
  if (param && ! param->is_string()) throw (ExceptionParametersBadType());
  return (param != NULL) ? param->value_string_ : deflt;
}

//----------------------------------------------------------------------

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
/// @param   index     Index into the list
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

//----------------------------------------------------------------------

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
/// @param   index     Index into the list
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

//----------------------------------------------------------------------

void Parameters::monitor_log (std::string parameter) throw()
{
  Param * param = parameter_(parameter);
  char buffer[80];
  sprintf (buffer,"Parameter %s:%s:%s = %s",
	   current_group_.c_str(),
	   current_subgroup_.c_str(),
	   parameter.c_str(),
	   param ? param->value_to_string().c_str() : "undefined");
  monitor_->print(buffer);
}
//----------------------------------------------------------------------

int Parameters::group_count() throw ()
{
  return parameter_tree_->size();
}

//----------------------------------------------------------------------

std::string Parameters::group(int group_index) throw ()
{
  return parameter_tree_->subgroup(group_index);
}

//----------------------------------------------------------------------

int Parameters::subgroup_count() throw ()
{
  if (parameter_tree_->subnode(current_group_) != 0) {
    return parameter_tree_->subnode(current_group_)->size();
  }
  return 0;
}

//----------------------------------------------------------------------

std::string Parameters::subgroup(int group_index) throw ()
{
  if (parameter_tree_->subnode(current_group_) != 0) {
    return parameter_tree_->subnode(current_group_)->subgroup(group_index);
  }
  return "";
}


//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

Param * Parameters::list_element_ (std::string parameter, int index) throw()
{
  Param * list = parameter_(parameter);
  Param * param = NULL;
  if (list == NULL) {
    char message [ ERROR_MESSAGE_LENGTH ];
    sprintf (message, 
	     "uninitialized parameter %s[%d] accessed\n",
	     parameter.c_str(),index);
    WARNING_MESSAGE("Parameters::list_element_",message);
  } else {
    int list_length = list->value_list_->size();
    if (list != NULL && 0 <= index && index < list_length ) {
      param =  (*(list->value_list_))[index];
    }
  }
  return param;
}

//----------------------------------------------------------------------

void Parameters::new_param_
(
 std::string group,
 std::string subgroup,
 std::string parameter,
 Param * param
 ) throw()
{

  std::string full_parameter = 
    group + ":" + subgroup + ":" + parameter;

  // Insert parameter into the parameter mapping
  // "Group:subgroup:parameter" -> "Value"

  parameter_map_     [full_parameter] = param;
    
  // Insert parameter into the parameter tree "Group" -> "subgroup"
  // -> "parameter"

  ParamNode * param_node = parameter_tree_;
  param_node = param_node->new_subnode(group);
  param_node = param_node->new_subnode(subgroup);
  param_node = param_node->new_subnode(parameter);
}
