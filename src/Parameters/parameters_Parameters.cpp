// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
/// @bug      Probable memory leaks
/// @todo     Add more info to Exception messages, e.g. parameter, expected type and actual type
/// @todo     Fix idiosynchratic ';' usage: currently must not be between top level groups, but must be between subgroups
/// @brief    Read in a parameter file and access parameter values

#include "cello.hpp"

#include "parameters.hpp"

// Parameters Parameters::instance_; // (singleton design pattern)

//----------------------------------------------------------------------

Parameters::Parameters(Monitor * monitor) 
  throw()
  : current_group_depth_(0),
    parameter_map_(),
    parameter_tree_(new ParamNode("Cello")),
    monitor_(monitor)
{
  if (! monitor_) monitor_ = Monitor::instance();

  for (int i=0; i<MAX_GROUP_DEPTH; i++) current_group_[i] = 0;
}

//----------------------------------------------------------------------

Parameters::Parameters(const char * file_name,
		       Monitor * monitor) 
  throw()
  : current_group_depth_(0),
    parameter_map_(),
    parameter_tree_(new ParamNode("Cello")),
    monitor_(monitor)
{
  if (! monitor_) monitor_ = Monitor::instance();
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

  monitor_->print ("[Parameter] read in %s",file_name);
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

  // "Previous" groups are empty
  int n_prev = 0;
  std::string group_prev[MAX_GROUP_DEPTH];

  // Initialize "current" groups as empty
  int n_curr = 0;
  std::string group_curr[MAX_GROUP_DEPTH];

  // Initialize indentation variables
  const std::string indent_amount = "    ";
  int indent_size = 4;
  std::string indent_string = "   ";
  int group_depth = 0;

  // Loop over parameters
  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {

    // Does the current parameter have a value?
    if (it_param->second) {

      // Determine groups of the current parameter
      int n_curr = extract_groups_(it_param->first,group_curr);

      // Determine the first group that differs, if any
      int i_group;
      for (i_group = 0; 
	   i_group < n_prev && i_group < n_curr &&
	     group_prev[i_group] == group_curr[i_group];
	   i_group++) {
	// (Intentionally blank)
      }

      // End old groups

      for (int i=i_group; i<n_prev; i++) {
	--group_depth;
	indent_string = indent_string.substr(indent_size,std::string::npos);
	fprintf (file_pointer, "%s}%c\n",indent_string.c_str(),
		 (group_depth==0) ? ' ' : ';' );
      }

      // Begin new groups

      for (int i=i_group; i<n_curr; i++) {
	fprintf (file_pointer,"%s%s {\n",indent_string.c_str(),
		 group_curr[i].c_str());
	indent_string = indent_string + indent_amount;
	++group_depth;
      }

      // Print parameter
      fprintf (file_pointer,"%s",indent_string.c_str());
      it_param->second->write(file_pointer,it_param->first);

      // Copy current groups to previous groups (inefficient)
      n_prev = n_curr;
      for (int i=0; i<n_prev; i++) {
	group_prev[i] = group_curr[i];
      }

    } else {
      char message [ ERROR_LENGTH ];
      sprintf (message, 
	       "uninitialized parameter %s accessed",
	       it_param->first.c_str());
      WARNING("Parameters::write",message);
    }
  }

  // End old groups

  for (int i=0; i<n_prev; i++) {
    indent_string = indent_string.substr(indent_size,std::string::npos);
    fprintf (file_pointer, "%s}\n",indent_string.c_str());
    --group_depth;
  }

  if (file_pointer != stdout) fclose(file_pointer);

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

double Parameters::value_float
( std::string parameter,
  double      deflt ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return floating point (double) parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_float()) throw ExceptionParametersBadType();
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  sprintf (deflt_string,"%#.15g",deflt);
  monitor_access_(parameter,deflt_string);
  return (param != NULL) ? param->get_float() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_float
( std::string parameter,
  double      value ) throw(ExceptionParametersBadType)
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);
  if (param) {
    if (! param->is_float()) throw ExceptionParametersBadType();
  } else {
    param = new Param;
    new_param_ (current_group_,parameter,param);
  }
  param->set_float_(value);
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

void Parameters::evaluate_float 
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
/// @param   result    Output array of evaluated floating point parameters values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{
  Param * param = parameter_(parameter);
  if (param && ! param->is_float_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_float(param->value_expr_,n,result,x,y,z,t);
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

double Parameters::list_value_float 
( int index,
  std::string parameter,
  double      deflt ) throw(ExceptionParametersBadType)
/// @param   index     Index of the floating point (double) list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return floating point (double) list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  if (param && ! param->is_float()) throw ExceptionParametersBadType();
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  sprintf (deflt_string,"%#.15g",deflt);
  monitor_access_(parameter,deflt_string,index);
  return (param != NULL) ? param->value_float_ : deflt;
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

void Parameters::list_evaluate_float 
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
/// @param   result    Output array of evaluated floating point expression list parameter element values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         Array of T values
{

  Param * param = list_element_(parameter,index);
  if (param && ! param->is_float_expr()) throw ExceptionParametersBadType();
  if (param != NULL) {
    param->evaluate_float(param->value_expr_,n,result,x,y,z,t);
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
  // Find the parameter node for the current list of groups
  ParamNode * param_node = parameter_tree_;
  for (int i=0; i<current_group_depth_; i++) {
    if (param_node->subnode(current_group_[i]) != 0) {
      param_node = param_node->subnode(current_group_[i]);
    }
  }
  int value = (param_node) ? param_node->size() : 0;
  return value;
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

void Parameters::group_set(int index, std::string group) throw()
{
  if (index >= MAX_GROUP_DEPTH) {
    char message [ ERROR_LENGTH ];
    sprintf (message, "group_set(%d,%s) index %d exceeds MAX_GROUP_DEPTH = %d",
	     index,group.c_str(),index,MAX_GROUP_DEPTH);
    ERROR("Parameters::group_set",message);
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

void Parameters::group_clear() throw ()
{
  for (int i=0; i<current_group_depth_; i++) {
    if (current_group_[i]) {
      free (current_group_[i]);
      current_group_[i] = NULL;
    }
  }
  current_group_depth_ = 0;
}
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

  std::string value;

  if ( param != NULL ) {
    value = param->value_to_string().c_str();
  } else {
    value = deflt_string + " [default]";
  }

  char index_string [MAX_PARAMETER_FILE_WIDTH] = "";

  if (index != -1) {
    sprintf (index_string,"[%d]",index);
  }

  monitor_->print("[Parameter] accessed %s%s %s",
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

  monitor_->print(buffer);
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

size_t Parameters::extract_groups_
(
 const std::string parameter, 
 std::string * group
)
{
  std::string p = parameter;
  size_t i_group=0;  // group index for group[]
  int i_stop  = p.find(":");
  while (i_stop != std::string::npos &&
	 i_group<MAX_GROUP_DEPTH) {
    group[i_group++] = p.substr(0,i_stop);
    p = p.substr(i_stop+1,std::string::npos);
    i_stop  = p.find(":");
  }
  return i_group;
}
//----------------------------------------------------------------------

void Parameters::new_param_
(
 char * group[],
 std::string parameter,
 Param * param
 ) throw()
{

  // @@@ Move into function: group[i] -> "group[0]:group[1]..."
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
