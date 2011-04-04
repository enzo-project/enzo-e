// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
/// @bug      Probable memory leaks
/// @todo     Add more info to Exception messages, e.g. parameter, expected type and actual type
/// @brief    Read in a parameter file and access parameter values

#include "cello.hpp"

#include "parameters.hpp"

// Parameters Parameters::instance_; // (singleton design pattern)

//----------------------------------------------------------------------

Parameters::Parameters() 
  throw()
  : current_group_depth_(0),
    parameter_map_(),
    parameter_tree_(new ParamNode("Parameters"))
{
  for (int i=0; i<MAX_GROUP_DEPTH; i++) current_group_[i] = 0;
}

//----------------------------------------------------------------------

Parameters::Parameters(const char * file_name ) 
  throw()
  : current_group_depth_(0),
    parameter_map_(),
    parameter_tree_(new ParamNode("Parameters"))
{
  for (int i=0; i<MAX_GROUP_DEPTH; i++) current_group_[i] = 0;
  read(file_name);
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

void Parameters::read ( const char * file_name )
/// @param    file_pointer An opened input parameter file or stdin
{

  FILE * file_pointer = fopen(file_name,"r");

  if ( file_pointer == NULL ) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Error opening parameter file '%s' for reading",file_name);
    ERROR("Parameters::read",buffer);
  }
  
  struct param_struct * parameter_list = cello_parameters_read(file_pointer);

  struct param_struct * node = parameter_list -> next; // skip sentinel
  struct param_struct * prev = node;

  while (node->type != enum_parameter_sentinel) {

    Param * param = new Param;

    param->set(node);

    new_param_(node->group,node->parameter,param);

    node = node->next;
    
    // free not delete since allocated in parse.y
    for (int i=0; prev->group[i] && i < MAX_GROUP_DEPTH; i++) {
      free (prev->group[i]);
    }
    free (prev);

    prev = node;
    
  }

  // assert: node->type == enum_parameter_sentinel

  free (node);
  node = NULL;

  fclose(file_pointer);

  Monitor::instance()->print ("[Parameter] read in %s",file_name);
}

//----------------------------------------------------------------------

void Parameters::write ( const char * file_name )
/// @param    file_pointer An opened output parameter file or stdout
{

  FILE * file_pointer = fopen(file_name,"w");

  if ( file_pointer == NULL ) {
    char buffer[ERROR_LENGTH];
    sprintf (buffer,"Error opening parameter file '%s' for writing",file_name);
    ERROR("Parameters::writing",buffer);
  }
  std::map<std::string,Param *>::iterator it_param;

  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {

    if (it_param->second) {
      it_param->second->write(file_pointer, it_param->first);
    } else {
      char message [ ERROR_LENGTH ];
      sprintf (message, 
	       "uninitialized parameter %s accessed",
	       it_param->first.c_str());
      WARNING("Parameters::write",message);
    }
  }

  fclose(file_pointer);

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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%d",deflt);
  monitor_access_(parameter,deflt_string);
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
    new_param_ (current_group_,parameter,param);
  }
  param->set_integer_(value);
  monitor_write_(parameter);
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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%g",deflt);
  monitor_access_(parameter,deflt_string);
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
    new_param_ (current_group_,parameter,param);
  }
  param->set_scalar_(value);
  monitor_write_(parameter);
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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%s",deflt ? "true" : "false");
  monitor_access_(parameter,deflt_string);
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
    new_param_ (current_group_,parameter,param);
  }
  param->set_logical_(value);
  monitor_write_(parameter);
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
  monitor_access_(parameter,deflt);
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
    new_param_ (current_group_,parameter,param);
  }
  param->set_string_(strdup(value));
  monitor_write_(parameter);
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
  // char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // sprintf_expression (param->value_expr_,deflt_string);
  // monitor_access_(parameter,deflt_string);
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
  // char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // sprintf_expression (param->value_expr_,deflt_string);
  // monitor_access_(parameter,deflt_string);
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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%d",deflt);
  monitor_access_(parameter,deflt_string,index);
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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%g",deflt);
  monitor_access_(parameter,deflt_string,index);
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
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%s",deflt ? "true" : "false");
  monitor_access_(parameter,deflt_string,index);
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
  monitor_access_(parameter,deflt,index);
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
  // char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // sprintf_expression (param->value_expr_,deflt_string);
  // monitor_access_(parameter,deflt_string,index);
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
  // char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // sprintf_expression (param->value_expr_,deflt_string);
  // monitor_access_(parameter,deflt_string,index);
}

//----------------------------------------------------------------------

std::string Parameters::group(int i) const throw()
{
  return (i < current_group_depth_) ? current_group_[i] : "";
}

//----------------------------------------------------------------------

int Parameters::group_depth() const throw()
{
  return current_group_depth_;
}

//----------------------------------------------------------------------

int Parameters::group_count() const throw()
{
  ParamNode * param_node = parameter_tree_;
  for (int i=0; i<current_group_depth_; i++) {
    if (param_node->subnode(current_group_[i]) != 0) {
      param_node = param_node->subnode(current_group_[i]);
    }
  }
  return (param_node) ? param_node->size() : 0;
}

//----------------------------------------------------------------------

void Parameters::group_push(std::string str) throw()
{
  if (current_group_depth_ < MAX_GROUP_DEPTH - 1) {
    current_group_[current_group_depth_++] = strdup(str.c_str());
  } else {
    char message [ ERROR_LENGTH ];
    sprintf (message, 
	     "Parameter grouping depth %d exceeds MAX_GROUP_DEPTH = %d",
	     current_group_depth_,MAX_GROUP_DEPTH);
    ERROR("Parameters::group_push",message);
  }
}

//----------------------------------------------------------------------

void Parameters::group_pop(std::string group) throw()
{
  if (current_group_depth_ > 0) {
    if (group != "" && group != current_group_[current_group_depth_-1]) {
      char message [ ERROR_LENGTH ];
      sprintf (message, "group_pop(%s) does not match group_push(%s)",
	       group.c_str(),current_group_[current_group_depth_-1]);
      WARNING("Parameters::group_pop",message);
    }
    --current_group_depth_;
    free (current_group_[current_group_depth_]);
    current_group_[current_group_depth_] = NULL;
  } else {
    ERROR("Parameters::group_pop",
	  "More calls to group_pop() than group_push()");
  }
}

//----------------------------------------------------------------------

void Parameters::set_group(int index, std::string group) throw()
{
  if (index >= MAX_GROUP_DEPTH) {
    char message [ ERROR_LENGTH ];
    sprintf (message, "set_group(%d,%s) index %d exceeds MAX_GROUP_DEPTH = %d",
	     index,group.c_str(),index,MAX_GROUP_DEPTH);
    ERROR("Parameters::set_group",message);
  }
  for (int i=index; i<current_group_depth_; i++) {
    if (current_group_[i]) {
      free (current_group_[i]);
      current_group_[i] = NULL;
    }
  }
  current_group_depth_ = index + 1;
  current_group_[index] = strdup(group.c_str());
}

//----------------------------------------------------------------------

// int Parameters::group_count() throw ()
// {
//   return parameter_tree_->size();
// }

// //----------------------------------------------------------------------

// std::string Parameters::group(int group_index) throw ()
// {
//   return parameter_tree_->subgroup(group_index);
// }

// //----------------------------------------------------------------------

// int Parameters::subgroup_count() throw ()
// {
//   if (parameter_tree_->subnode(current_group_) != 0) {
//     return parameter_tree_->subnode(current_group_)->size();
//   }
//   return 0;
// }

// //----------------------------------------------------------------------

// std::string Parameters::subgroup(int group_index) throw ()
// {
//   if (parameter_tree_->subnode(current_group_) != 0) {
//     return parameter_tree_->subnode(current_group_)->subgroup(group_index);
//   }
//   return "";
// }

//----------------------------------------------------------------------

parameter_enum Parameters::type
( std::string  parameter) throw()
/// @param   parameter Parameter name
/// @return  Return type of the given parameter
{
  Param * param = parameter_(parameter);
  return param ? param->type() : parameter_unknown ;
}

//----------------------------------------------------------------------

parameter_enum Parameters::list_type
( 
 int index,
 std::string  parameter
 ) throw()
/// @param   index Index into the list
/// @param   parameter Parameter name
/// @return  Return type of the given parameter
{
  Param * list = parameter_(parameter);
  Param * param = NULL;
  if (list != NULL) {
    int list_length = list->value_list_->size();
    if (list != NULL && 0 <= index && index < list_length ) {
      param =  (*(list->value_list_))[index];
    }
  }
  return param ? param->type() : parameter_unknown ;
}

//======================================================================

void Parameters::monitor_access_ 
(
 std::string parameter,
 std::string deflt_string,
 int index
) throw()
{
  Param * param = 0;

  if (index == -1) {
    // not a list element
    param = parameter_(parameter);
  } else {
    // a list element
    param = list_element_(parameter,index);
  }
  std::string value = param ? 
    param->value_to_string().c_str() : 
    std::string("[" + deflt_string + "]");

  char index_string [MAX_PARAMETER_FILE_WIDTH] = "";

  if (index != -1) {
    sprintf (index_string,"[%d]",index);
  }

  Monitor * monitor = Monitor::instance();
  monitor->print("[Parameter] accessed %s%s %s",
		 parameter_name_(parameter).c_str(),
		 index_string,
		 value.c_str());
}

//----------------------------------------------------------------------

void Parameters::monitor_write_ (std::string parameter) throw()
{
  Param * param = parameter_(parameter);
  char buffer[MONITOR_LENGTH];
  sprintf (buffer,"Parameter write %s = %s",
	   parameter_name_(parameter).c_str(),
	   param ? param->value_to_string().c_str() : "[undefined]");

  Monitor * monitor = Monitor::instance();
  monitor->print(buffer);
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
  if (list != NULL) {
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
 char * group[],
 std::string parameter,
 Param * param
 ) throw()
{

  std::string full_parameter = "";
  for (int i=0; group[i] != 0 && i < MAX_GROUP_DEPTH; i++) {
    full_parameter = full_parameter + group[i] + ":";
  }
  full_parameter = full_parameter + parameter;

  // Insert parameter into the parameter mapping
  // "Group:group[0]:group[1]:...:parameter" -> "Value"

  parameter_map_     [full_parameter] = param;
    
  // Insert parameter into the parameter tree "Group" -> "group[0]" ->
  // "group[1]" -> ... -> "parameter"

  ParamNode * param_node = parameter_tree_;
  for (int i=0; group[i] != NULL && i < MAX_GROUP_DEPTH; i++) {
    param_node = param_node->new_subnode(group[i]);
  }
  param_node = param_node->new_subnode(parameter);
}
