// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:38:43 PDT 2009
/// @brief    Read in a parameter file and access parameter values

#include "cello.hpp"

#include "parameters.hpp"

Parameters g_parameters;

// #define TEMP_BYPASS_DESTRUCTOR

//----------------------------------------------------------------------

Parameters::Parameters(Monitor * monitor) 
  throw()
  : current_group_(),
    parameter_map_(),
    parameter_tree_(new ParamNode("Cello")),
    monitor_(monitor),
    lmonitor_(true)
{
  if (! monitor_) lmonitor_ = false;
}

//----------------------------------------------------------------------

Parameters::Parameters(const char * file_name,
		       Monitor * monitor) 
  throw()
  : current_group_(),
    parameter_map_(),
    parameter_tree_(new ParamNode("Cello")),
    monitor_(monitor),
    lmonitor_(true)
{
  if (! monitor_) lmonitor_ = false;
  read(file_name);
}

//----------------------------------------------------------------------

Parameters::~Parameters()
///
{
  // Iterate over all parameters, deleting their values


#ifdef TEMP_BYPASS_DESTRUCTOR
  for (auto it_param =  parameter_map_.begin();
       it_param != parameter_map_.end();
       ++it_param) {
    delete it_param->second;
  }
#endif  
  delete parameter_tree_;
}

//----------------------------------------------------------------------

void Parameters::pup (PUP::er &p)
{
  TRACEPUP;

  p | current_group_;

  int n = 0;
  if (!p.isUnpacking()) {
    // Figure out size
    for (auto it_param =  parameter_map_.begin();
	 it_param != parameter_map_.end();
	 ++it_param) {
      n++;
    }
  }
  p | n;
  if (!p.isUnpacking()) {
    for (auto it_param =  parameter_map_.begin();
	 it_param != parameter_map_.end();
	 ++it_param) {
      std::string name = it_param->first;
      Param * param = it_param->second;
      p | name;
      int lparam = (param != NULL);
      p | lparam;
      ASSERT("Parameters::pup",
             "param is expected to be non-null",
             ((! lparam) && (param == NULL)) ||
             (  lparam  && (param != NULL)));
      if (lparam) p | *param;
    } 
  } else {
    for (int i=0; i<n; i++) {
      std::string name;
      p | name;
      int lparam=0;
      p | lparam;
      Param * param = NULL;
      if (lparam) {
	param = new Param;
	p | *param;
      }
      parameter_map_[name] = param;
    }
  }
}

//----------------------------------------------------------------------

void Parameters::read ( const char * file_name )
/// @param    file_name An opened input parameter file or stdin
{

  FILE * file_pointer = fopen(file_name,"r");

  ASSERT1("Parameters::read",
	  "Error opening parameter file '%s' for reading",
	  file_name,
	  ( file_pointer != NULL ));

  struct param_struct * parameter_list = 
    cello_parameters_read(file_name,file_pointer);

  struct param_struct * node = parameter_list -> next; // skip sentinel
  struct param_struct * prev = node;

  while (node->type != enum_parameter_sentinel) {

    Param * param = new Param;

    param->set(node);

    std::string full_parameter = "";
    for (int i=0;  (i < MAX_GROUP_DEPTH) && node->group[i] != 0; i++) {
      full_parameter = full_parameter + node->group[i] + ":";
    }

    full_parameter = full_parameter + node->parameter;

    new_param_(full_parameter,param);

    node = node->next;
    
    // free not delete since allocated in parse.y
    for (int i=0; (i < MAX_GROUP_DEPTH) && prev->group[i]; i++) {
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

void Parameters::write ( const char * file_name, int type )
/// @param  file_name   An opened output parameter file or stdout
{
  FILE * fp = fopen(file_name,"w");
  ASSERT1("Parameters::write",
	  "Error opening parameter file '%s' for writing",
	  file_name,
	  ( fp != NULL ) );

  write (fp,type);
  fclose(fp);
}

//----------------------------------------------------------------------

void Parameters::write ( FILE * fp, int type )
/// @param  file_name   An opened output parameter file or stdout
{
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

  for (auto it_param =  parameter_map_.begin();
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
	if (type != param_write_monitor) {
	  fprintf (fp, "%s}%c\n",indent_string.c_str(),
		   (group_depth==0) ? '\n' : ';' );
	}
      }

      // Begin new groups

      for (int i=i_group; i<n_curr; i++) {
	if (type != param_write_monitor) {
	  const char * format =
	    (type==param_write_libconfig) ? "%s%s : {\n" : "%s%s {\n";
	  fprintf (fp, format ,indent_string.c_str(),
		   group_curr[i].c_str());
	}
	indent_string = indent_string + indent_amount;
	++group_depth;
      }

      // Print parameter
      if (type == param_write_monitor) {
	Monitor monitor;
	monitor.set_mode(monitor_mode_all);
	// display Monitor prefix if full_names
	monitor.write(fp,"Parameters","");
	for (int i=0; i < n_curr; i++) {
	  fprintf (fp,"%s:",group_curr[i].c_str());
	}
      } else {	
	fprintf (fp,"%s",indent_string.c_str());
      }
      it_param->second->write(fp,it_param->first,type);

      // Copy current groups to previous groups (inefficient)
      n_prev = n_curr;
      for (int i=0; i<n_prev; i++) {
	group_prev[i] = group_curr[i];
      }

    }
  }

  // End old groups

  for (int i=0; i<n_prev; i++) {
    indent_string = indent_string.substr(indent_size,std::string::npos);
    if (type != param_write_monitor) {
      fprintf (fp, "%s}\n",indent_string.c_str());
    }
    --group_depth;
  }

}

//----------------------------------------------------------------------

int Parameters::value_integer 
( std::string parameter,
  int         deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return integer parameter value if it exists, deflt if not
{
  Param * param = this->param(parameter);

  ASSERT1 ("Parameters::value_integer",
	   "Parameter %s is not an integer", parameter.c_str(),
	   ( ! param || param->is_type(parameter_integer)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,"%d",deflt);
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
  Param * param = this->param(parameter_name_(parameter));

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
  Param * param = this->param(parameter);

  ASSERT1 ("Parameters::value_float",
	   "Parameter %s is not a float", parameter.c_str(),
	   ( ! param || param->is_type(parameter_float)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,FLOAT_FORMAT,deflt);
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
  Param * param = this->param(parameter);

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
  Param * param = this->param(parameter);

  ASSERT2 ("Parameters::value_logical",
	   "Parameter %s is type %d not a logical",
	   parameter.c_str(),param->type(),
	   ( ! param || param->is_type(parameter_logical)));

  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,
            "%s",deflt ? "true" : "false");
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
  Param * param = this->param(parameter);

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

std::string Parameters::value_string 
( std::string parameter,
  std::string deflt ) throw()
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string parameter value if it exists, deflt if not
{
  Param * param = this->param(parameter);

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
  Param * param = this->param(parameter);

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

Expression Parameters::value_Expression
( std::string parameter) throw()
/// @param   parameter Parameter name
/// @return  Return Expression parameter if it exists
{
  Param * param = this->param (parameter);
  ASSERT1 ("Parameters::value_string",
	   "Parameter %s[%d] does not exist",
	   parameter.c_str(), param!= nullptr);
  ASSERT1 ("Parameters::list_value_string",
	   "Parameter %s is not an expression",
	   parameter.c_str(), ( param->is_type(parameter_float_expr) ||
                                param->is_type(parameter_logical_expr) ));
  monitor_access_(parameter,"");
  param->set_accessed();
  return this->construct_or_rebuild_Expression_(parameter, -1);
}

//----------------------------------------------------------------------

int Parameters::list_length(std::string parameter)
/// @param   parameter Parameter name
{
  Param * param = this->param(parameter);
  ASSERT1 ("Parameters::list_length",
	   "Parameter %s is not a list", full_name(parameter).c_str(),
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
  Param * param = this->param(parameter,index);
  ASSERT2 ("Parameters::list_value_integer",
	   "Parameter %s[%d] is not an integer", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_integer)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,"%d",deflt);
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
  Param * param = this->param(parameter,index);
  ASSERT2 ("Parameters::list_value_float",
	   "Parameter %s[%d] is not a float", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_float)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  // '#' format character forces a decimal point
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,FLOAT_FORMAT,deflt);
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
  Param * param = this->param(parameter,index);
  ASSERT2 ("Parameters::list_value_logical",
	   "Parameter %s[%d] is not a logical", 
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_logical)));
  char deflt_string[MAX_PARAMETER_FILE_WIDTH];
  snprintf (deflt_string,MAX_PARAMETER_FILE_WIDTH,
            "%s",deflt ? "true" : "false");
  monitor_access_(parameter,deflt_string,index);
  return (param != NULL) ? param->value_logical_ : deflt;
}

//----------------------------------------------------------------------

std::string Parameters::list_value_string 
( int index,
  std::string parameter,
  std::string deflt ) throw()
/// @param   index     Index of the string list parameter element
/// @param   parameter Parameter name
/// @param   deflt     Default parameter value
/// @return  Return string list parameter element value if it exists, deflt if not
{
  Param * param = this->param (parameter,index);
  ASSERT2 ("Parameters::list_value_string",
	   "Parameter %s[%d] is not a string",
	   parameter.c_str(),index,
	   ( ! param || param->is_type(parameter_string)));
  monitor_access_(parameter,deflt,index);
  return (param != NULL) ? param->value_string_ : deflt;
}

//----------------------------------------------------------------------

Expression Parameters::list_value_Expression
( int index,
  std::string parameter) throw()
/// @param   index     Index of the Expression list parameter element
/// @param   parameter Parameter name
/// @return  Return Expression list parameter element value if it exists
{
  Param * param = this->param (parameter,index);
  ASSERT2 ("Parameters::list_value_string",
	   "Parameter %s[%d] does not exist",
	   parameter.c_str(),index, param!= nullptr);
  ASSERT2 ("Parameters::list_value_string",
	   "Parameter %s[%d] is not an expression",
	   parameter.c_str(), index, (param->is_type(parameter_float_expr) ||
				      param->is_type(parameter_logical_expr)));
  monitor_access_(parameter,"",index);
  param->set_accessed();
  return this->construct_or_rebuild_Expression_(parameter, index);
}

//----------------------------------------------------------------------

void Parameters::set_list_length 
(
 std::string parameter,
 int         length
 )
{
  Param * param = this->param(parameter);

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
  Param * param = this->param(parameter,index);

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

  Param * param = this->param(parameter,index);

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
  Param * param = this->param(parameter,index);

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
  Param * param = this->param(parameter,index);

  if ( ! param ) {
    param = new Param;
    new_param_ (parameter_name_(parameter),param);
  }

  param->set_string_(strdup(value));
}

//----------------------------------------------------------------------

std::string Parameters::group(int i) const throw()
{
  return (i < current_group_.size()) ? current_group_[i] : "";
}

//----------------------------------------------------------------------

int Parameters::group_depth() const throw()
{
  return current_group_.size();
}

//----------------------------------------------------------------------

int Parameters::group_count() const throw()
{
  // Find the parameter node for the current list of groups
  ParamNode * param_node = parameter_tree_;
  for (int i=0; i<current_group_.size(); i++) {
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
  int n = current_group_.size();
  current_group_.resize(n + 1);
  current_group_[n] = str;
}

//----------------------------------------------------------------------

void Parameters::group_pop(std::string group) throw()
{
  int n = current_group_.size();
  
  ASSERT("Parameters::group_pop",
	"More calls to group_pop() than group_push()",
	 (n > 0));

  if (group != "" && group != current_group_[n-1]) {
    WARNING2("Parameters::group_pop",
	     "group_pop(%s) does not match group_push(%s)",
	     group.c_str(),current_group_[n-1].c_str());
  }

  current_group_.resize(n - 1);
  
}

//----------------------------------------------------------------------

void Parameters::group_set(int index, std::string group) throw()
{
  current_group_.resize(index+1);
  current_group_[index] = group;
}

//----------------------------------------------------------------------

void Parameters::group_clear() throw ()
{
  current_group_.resize(0);
}

//----------------------------------------------------------------------

std::vector<std::string> Parameters::leaf_parameter_names() const throw()
{
  std::vector<std::string> out;
  std::string prefix = full_name("");
  std::size_t prefix_size = prefix.size();
  
  for (auto it_param =  parameter_map_.lower_bound(prefix);
       it_param != parameter_map_.end();
       ++it_param) {

    // If the current key doesn't share the prefix, abort the search
    if (((it_param->first).compare(0, prefix_size, prefix) != 0) ||
	((it_param->first).size() <= prefix_size) ) {
      break;
    }
    std::string suffix = (it_param->first).substr(prefix_size,
						  std::string::npos);
    if (suffix.find(':') != std::string::npos){
      // this isn't a leaf parameter
      continue;
    } else {
      out.push_back(suffix);
    }
  }
  return out;
}

//----------------------------------------------------------------------

std::string Parameters::full_name(std::string parameter) const throw()
{
  int n = current_group_.size();
  std::string full_name = "";
  for (int i=0; i<n; i++) {
    full_name = full_name + current_group_[i] + ":";
  }
  full_name = full_name + parameter;
  return full_name;

}
//----------------------------------------------------------------------

parameter_type Parameters::type
( std::string  parameter) throw()
/// @param   parameter Parameter name
/// @return  Return type of the given parameter
{
  Param * param = this->param(parameter);
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
  Param * list = this->param(parameter);
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
    param = this->param(parameter);
  } else {
    // a list element
    param = this->param(parameter,index);
  }

  std::string value;

  if ( param != NULL ) {
    value = param->value_to_string(param_write_monitor);
  } else {
    value = deflt_string + " [default]";
  }

  std::string index_string = 
    std::string("[") + std::to_string(index) + std::string("]");

  std::string buffer = std::string("accessed ") +
    parameter_name_(parameter) + index_string + " " + value;

  if (lmonitor_) monitor_->print_verbatim("Parameters",buffer.data());
}

//----------------------------------------------------------------------

void Parameters::monitor_write_ (std::string parameter) throw()
{
  Param * param = this->param(parameter);
  char buffer[MONITOR_LENGTH];
  snprintf (buffer,MONITOR_LENGTH,"Parameter write %s = %s",
	   parameter_name_(parameter).c_str(),
	   param ?
	   param->value_to_string(param_write_monitor).c_str() : "[undefined]");

  if (lmonitor_) monitor_->print_verbatim("Parameters",buffer);
}

//----------------------------------------------------------------------

int Parameters::readline_ 
( FILE*  fp, 
  char * buffer, 
  int    buffer_length ) 
  throw()
/// @param   fp                 the opened input file
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
    c = fgetc(fp);
    buffer[i] = (std::numeric_limits<char>::min() <= c && 
		 c <= std::numeric_limits<char>::max()) ? c : '\0';
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

/// Return the const Param pointer for the specified parameter
std::pair<const Param*, const Param*> Parameters::const_param_
(std::string parameter, int index) const
{
  const Param *list_ptr, *param_ptr;
  auto search = parameter_map_.find(parameter_name_(parameter));
  if (search == parameter_map_.end()){
    list_ptr = nullptr;
    param_ptr = nullptr;
  } else if (index == -1){
    list_ptr = nullptr;
    param_ptr = search->second;
  } else {
    list_ptr = search->second;
    param_ptr = nullptr;
    int list_length = list_ptr->value_list_->size();
    if (0 <= index && index < list_length ) {
      param_ptr = (*(list_ptr->value_list_))[index];
    }
  }
  return std::make_pair(param_ptr, list_ptr);
}

//----------------------------------------------------------------------

/// Return the Param pointer for the specified parameter
Param * Parameters::param (std::string parameter, int index)
{
  std::pair<const Param*, const Param*> ptr_pair = const_param_(parameter,
								index);

  // Casting const Param* to Param * is okay as long as parameter_map_ is
  // defined as a non-const instance of std::map<std::string, Param *>
  Param* param_ptr = const_cast<Param*>(ptr_pair.first);
  Param* list_ptr = const_cast<Param*>(ptr_pair.second);

  if (list_ptr != nullptr){ list_ptr->set_accessed(); }
  return param_ptr;
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
  for (auto it_param =  parameter_map_.begin();
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

Expression Parameters::construct_or_rebuild_Expression_
(const std::string & param_name, int param_index) const throw()
{
  std::string full_name = this->parameter_name_(param_name);

  std::pair<const Param*, const Param*> ppair = const_param_(full_name,
							     param_index);
  const Param* param_ptr = ppair.first;
  // ppair.second is just the list holding ppair.first (if param_index != -1)

  // sanity check (this is somewhat redundant on construction)
  ASSERT("Parameters::construct_or_rebuild_Expression_",
	 "The parameter is not an expression",
	 ( (param_ptr != nullptr) &&
	   (param_ptr->is_type(parameter_float_expr) ||
	    param_ptr->is_type(parameter_logical_expr)) ) );
  const struct node_expr * value_expr = param_ptr->value_expr_;

  return Expression(value_expr, param_ptr->is_type(parameter_float_expr),
		    full_name, param_index);
}

//----------------------------------------------------------------------
