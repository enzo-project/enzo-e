// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
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
  if (! monitor_) lmonitor_ = false;

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
  if (! monitor_) lmonitor_ = false;
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
  for (int i=0; i<current_group_depth_; i++) {
    free (current_group_[i]);
    current_group_[i] = 0;
  }
}

//----------------------------------------------------------------------

void Parameters::read ( const char * file_name )
/// @param    file_pointer An opened input parameter file or stdin
{

  FILE * file_pointer = fopen(file_name,"r");

  if ( file_pointer == NULL ) {
    ERROR1("Parameters::read",
	   "Error opening parameter file '%s' for reading",
	   file_name);
  }
  
  struct param_struct * parameter_list = 
    cello_parameters_read(file_name,file_pointer);

  struct param_struct * node = parameter_list -> next; // skip sentinel
  struct param_struct * prev = node;

  while (node->type != enum_parameter_sentinel) {

    Param * param = new Param;

    param->set(node);

    std::string full_parameter = "";
    for (int i=0; node->group[i] != 0 && i < MAX_GROUP_DEPTH; i++) {
      full_parameter = full_parameter + node->group[i] + ":";
    }
    full_parameter = full_parameter + node->parameter;

    new_param_(full_parameter,param);

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

  if (lmonitor_) monitor_->print ("Parameters","read in %s",file_name);
}

//----------------------------------------------------------------------

void Parameters::write ( const char * file_name )
/// @param    file_pointer An opened output parameter file or stdout
{

  FILE * file_pointer = fopen(file_name,"w");

  if ( file_pointer == NULL ) {
    ERROR1("Parameters::write",
	   "Error opening parameter file '%s' for writing",
	   file_name);
  }

  // "Previous" groups are empty
  int n_prev = 0;
  std::string group_prev[MAX_GROUP_DEPTH];

  // Initialize "current" groups as empty
  std::string group_curr[MAX_GROUP_DEPTH];

  // Initialize indentation variables
  const std::string indent_amount = "    ";
  int indent_size = 4;
  std::string indent_string = " ";
  int group_depth = 0;

  // Loop over parameters

  std::map<std::string,Param *>::iterator it_param;

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
	indent_string = indent_string.substr(0, indent_string.size()-indent_size);
	fprintf (file_pointer, "%s}%c\n",indent_string.c_str(),
		 (group_depth==0) ? '\n' : ';' );
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
      // OK--just means parameter default was used
      // WARNING1("Parameters::write",
      // 	       "uninitialized parameter %s accessed",
      // 	       it_param->first.c_str());
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
  int         deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::value_integer",
	   "Parameter %s is not an integer", parameter.c_str(),
	   ( ! param || param->is_type(parameter_integer)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%d",deflt);
  monitor_access_(parameter,deflt_string);
  return (param != NULL) ? param->get_integer() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_integer 
( std::string parameter,
  int         value ) throw()
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter_name_(parameter));

  ASSERT1 ("Parameters::set_integer",
	   "Parameter %s is not an integer", parameter.c_str(),
	   ( ! param || param->is_type(parameter_integer)));

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_integer_(value);
  monitor_write_(parameter);
}

//----------------------------------------------------------------------

double Parameters::value_float
( std::string parameter,
  double      deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return floating point (double) parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::value_float",
	   "Parameter %s is not a float", parameter.c_str(),
	   ( ! param || param->is_type(parameter_float)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  sprintf (deflt_string,FLOAT_FORMAT,deflt);
  monitor_access_(parameter,deflt_string);
  return (param != NULL) ? param->get_float() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_float
( std::string parameter,
  double      value ) throw()
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::set_float",
	   "Parameter %s is not a float", parameter.c_str(),
	   (!param || param->is_type(parameter_float)));

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_float_(value);
  monitor_write_(parameter);
}

//----------------------------------------------------------------------

bool Parameters::value_logical 
( std::string parameter,
  bool        deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return logical parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::value_logical",
	   "Parameter %s is not a logical", parameter.c_str(),
	   ( ! param || param->is_type(parameter_logical)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%s",deflt ? "true" : "false");
  monitor_access_(parameter,deflt_string);
  return (param != NULL) ? param->get_logical() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_logical
( std::string parameter,
  bool        value ) throw()
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::set_logical",
	   "Parameter %s is not logical", parameter.c_str(),
	   (!param || param->is_type(parameter_logical)));

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_logical_(value);
  monitor_write_(parameter);
}

//----------------------------------------------------------------------

const char * Parameters::value_string 
( std::string  parameter,
  const char * deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string parameter value if it exists, deflt if not
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::value_string",
	   "Parameter %s is not a string", parameter.c_str(),
	   ( ! param || param->is_type(parameter_string)));

  monitor_access_(parameter,deflt);
  return (param != NULL) ? param->get_string() : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_string
( std::string  parameter,
  const char * value ) throw()
/// @param   parameter Parameter name
/// @param   value     Value to set the parameter
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::set_string_string",
	   "Parameter %s is not a string", parameter.c_str(),
	   ( ! param || param->is_type(parameter_string)));

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
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
 double    t) throw()
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated floating point parameters values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of x values
/// @param   y         Array of y values
/// @param   z         Array of z values
/// @param   t         t value
{
  Param * param = parameter_(parameter);
  ASSERT1 ("Parameters::evaluate_float",
	   "Parameter %s is not a floating-point expression", parameter.c_str(),
	   ( ! param || param->is_type(parameter_float_expr)));
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
 double    t) throw()
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated logical parameters values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         T value
{
  Param * param = parameter_(parameter);
  ASSERT1 ("Parameters::evaluate_logical",
	   "Parameter %s is not a logical expression", parameter.c_str(),
	   (! param || param->is_type(parameter_logical_expr)));
  if (param != NULL) {
    param->evaluate_logical(param->value_expr_,n,result,x,y,z,t);
  } else {
    WARNING("Parameters::evaluate_logical",
	    "param is NULL but deflt not set");
  //   for (int i=0; i<n; i++) result[i] = deflt[i];
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
  ASSERT1 ("Parameters::list_length",
	   "Parameter %s is not a list", parameter.c_str(),
	   ( ! param || param->is_type(parameter_list)));
  return (param != NULL) ? (param->value_list_)->size() : 0;
}

//----------------------------------------------------------------------

int Parameters::list_value_integer 
( int index,
  std::string parameter,
  int         deflt ) throw()
/// @param   index     Index of the integer list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  ASSERT2 ("Parameters::list_value_integer",
	   "Parameter %s[%d] is not an integer", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_integer)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%d",deflt);
  monitor_access_(parameter,deflt_string,index);
  return (param != NULL) ? param->value_integer_ : deflt;
}

//----------------------------------------------------------------------

double Parameters::list_value_float 
( int index,
  std::string parameter,
  double      deflt ) throw()
/// @param   index     Index of the floating point (double) list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return floating point (double) list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  ASSERT2 ("Parameters::list_value_float",
	   "Parameter %s[%d] is not a float", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_float)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  sprintf (deflt_string,FLOAT_FORMAT,deflt);
  monitor_access_(parameter,deflt_string,index);
  return (param != NULL) ? param->value_float_ : deflt;
}

//----------------------------------------------------------------------

bool Parameters::list_value_logical 
( int index,
  std::string parameter,
  bool        deflt ) throw()
/// @param   index     Index of the logical list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return logical list parameter element value if it exists, deflt if not
{
  Param * param = list_element_(parameter,index);
  ASSERT2 ("Parameters::list_value_logical",
	   "Parameter %s[%d] is not a logical", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_logical)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  sprintf (deflt_string,"%s",deflt ? "true" : "false");
  monitor_access_(parameter,deflt_string,index);
  return (param != NULL) ? param->value_logical_ : deflt;
}

//----------------------------------------------------------------------

const char * Parameters::list_value_string 
( int index,
  std::string   parameter,
  const char *  deflt ) throw()
/// @param   index     Index of the string list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string list parameter element value if it exists, deflt if not
{
  Param * param = list_element_ (parameter,index);
  ASSERT2 ("Parameters::list_value_string",
	   "Parameter %s[%d] is not a string",
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_string)));
  monitor_access_(parameter,deflt,index);
  return (param != NULL) ? param->value_string_ : deflt;
}

//----------------------------------------------------------------------

void Parameters::set_list_length 
(
 std::string parameter,
 int         length
 )
{
  Param * param = parameter_(parameter);

  ASSERT1 ("Parameters::set_list_length",
	   "Parameter %s is not a list", parameter.c_str(),
	   ( ! param || param->is_type(parameter_list)));

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
    param->type_ = parameter_list;
    typedef std::vector<class Param *> list_type;
    param->value_list_ = new list_type;
  }
  param->value_list_->resize(length,0);
  for (int i =0; i<length; i++) {
    if ((*(param->value_list_))[i]==0)
      (*(param->value_list_))[i]=new Param;
  }
}

//----------------------------------------------------------------------

void Parameters::set_list_integer 
(
 int         index,
 std::string parameter,
 int         value
) throw()
{
  Param * param = list_element_(parameter,index);

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_integer_(value);
}

//----------------------------------------------------------------------


void Parameters::set_list_float
(
 int         index,
 std::string parameter,
 double      value
) throw()
{

  Param * param = list_element_(parameter,index);

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_float_(value);
}

//----------------------------------------------------------------------


void Parameters::set_list_logical
(
 int         index,
 std::string parameter,
 bool        value
) throw()
{
  Param * param = list_element_(parameter,index);

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_logical_(value);
}

//----------------------------------------------------------------------


void Parameters::set_list_string
(
 int          index,
 std::string  parameter,
 const char * value
) throw()
{
  Param * param = list_element_(parameter,index);

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_string_(strdup(value));
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
 double    t
 ) throw()
/// @param   index     Index into the list
/// @param   parameter Parameter name
/// @param   n         Length of variable arrays
/// @param   result    Output array of evaluated floating point expression list parameter element values if it exists, or deflt if not
/// @param   deflt     Array of default values
/// @param   x         Array of X values
/// @param   y         Array of Y values
/// @param   z         Array of Z values
/// @param   t         T value
{

  Param * param = list_element_(parameter,index);
  ASSERT2 ("Parameters::list_evaluate_float",
	   "Parameter %s[%d] is not a floating-point expression",
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_float_expr)));
  if (param != NULL) {
    param->evaluate_float(param->value_expr_,n,result,x,y,z,t);
  // } else {
  //   for (int i=0; i<n; i++) result[i] = deflt[i];
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
 double    t) throw()
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
  ASSERT2 ("Parameters::list_evaluate_logical",
	   "Parameter %s[%d] is not a logical",
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_logical_expr)));
  if (param != NULL) {
    param->evaluate_logical(param->value_expr_,n,result,x,y,z,t);
  } else {
    WARNING("Parameters::list_evaluate_logical",
	    "param is NULL but deflt not set");
  //   for (int i=0; i<n; i++) result[i] = deflt[i];
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
    ERROR2("Parameters::group_push",
	   "Parameter grouping depth %d exceeds MAX_GROUP_DEPTH = %d",
	   current_group_depth_,MAX_GROUP_DEPTH);
  }
}

//----------------------------------------------------------------------

void Parameters::group_pop(std::string group) throw()
{
  if (current_group_depth_ > 0) {
    if (group != "" && group != current_group_[current_group_depth_-1]) {
      WARNING2("Parameters::group_pop",
	       "group_pop(%s) does not match group_push(%s)",
	       group.c_str(),current_group_[current_group_depth_-1]);
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
    ERROR4("Parameters::group_set",
	   "group_set(%d,%s) index %d exceeds MAX_GROUP_DEPTH = %d",
	   index,group.c_str(),index,MAX_GROUP_DEPTH);
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

//----------------------------------------------------------------------

parameter_type Parameters::type
( std::string  parameter) throw()
/// @param   parameter Parameter name
/// @return  Return type of the given parameter
{
  Param * param = parameter_(parameter);
  return param ? param->type() : parameter_unknown ;
}

//----------------------------------------------------------------------

parameter_type Parameters::list_type
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
  if (! lmonitor_) return;

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
    value = param->value_to_string();
  } else {
    value = deflt_string + " [default]";
  }

  char index_string [MAX_PARAMETER_FILE_WIDTH] = "";

  if (index != -1) {
    sprintf (index_string,"[%d]",index);
  }

  char buffer[MAX_PARAMETER_FILE_WIDTH];
  sprintf (buffer,"accessed %s%s %s",
		  parameter_name_(parameter).c_str(),
		  index_string,
	   value.c_str());
  if (lmonitor_) monitor_->print_verbatim("Parameters",buffer);
}

//----------------------------------------------------------------------

void Parameters::monitor_write_ (std::string parameter) throw()
{
  Param * param = parameter_(parameter);
  char buffer[MONITOR_LENGTH];
  sprintf (buffer,"Parameter write %s = %s",
	   parameter_name_(parameter).c_str(),
	   param ? param->value_to_string().c_str() : "[undefined]");

  if (lmonitor_) monitor_->print_verbatim("Parameters",buffer);
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

  ASSERT ("Parameters::readline_",
	  "Input file line is too long",
	  (i+1 >= buffer_length-1));

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
    list->set_accessed();
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
  size_t i_stop  = p.find(":");
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
 std::string full_parameter,
 Param * param
 ) throw()
{
  parameter_map_ [full_parameter] = param;
    
  std::string groups[MAX_GROUP_DEPTH];

  int num_groups = extract_groups_(full_parameter,groups);

  ParamNode * param_node = parameter_tree_;
  for (int i=0; i<num_groups; i++) {
    param_node = param_node->new_subnode(groups[i]);
  }
  size_t ic = full_parameter.rfind(":");
  std::string parameter = full_parameter;
  if (ic != std::string::npos) {
    parameter = full_parameter.substr(ic+1,std::string::npos);
  }
  param_node = param_node->new_subnode(full_parameter);
}

//----------------------------------------------------------------------

void Parameters::check()
{
  std::map<std::string,Param *>::iterator it_param;

  for (it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {
    if (it_param->second && ! it_param->second->accessed()) {
      WARNING1 ("Parameters::check()",
		"Parmeter \"%s\" not accessed",
		it_param->first.c_str());
    }
  }
}

//----------------------------------------------------------------------
