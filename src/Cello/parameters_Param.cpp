// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Param.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sun Oct 11 15:02:08 PDT 2009
/// @brief    Implementation of the Param class

#include "cello.hpp"

#include "parameters.hpp"

const std::map<std::string, double (*) (double)> Param::function_map =
  {
    {"acos",acos},
    {"acosh", acosh},
    {"asin", asin},
    {"asinh", asinh},
    {"atan", atan},
    {"atanh", atanh},
    {"cbrt", cbrt},
    {"ceil", ceil},
    {"cos", cos},
    {"cosh", cosh},
    {"erfc", erfc},
    {"erf", erf},
    {"exp", exp},
    {"expm1", expm1},
    {"fabs", fabs},
    {"floor", floor},
    {"j0", j0},
    {"j1", j1},
    {"lgamma", lgamma},
    {"log10", log10},
    {"log1p", log1p},
    {"logb", logb},
    {"log", log},
    {"sin", sin},
    {"sinh", sinh},
    {"sqrt", sqrt},
    {"tan", tan},
    {"tanh", tanh},
    {"y0", y0},
    {"y1", y1},
    {"rint", rint},
  };

//----------------------------------------------------------------------

void Param::pup (PUP::er &p)
{
  TRACEPUP;
  const bool up = p.isUnpacking();
  // NOTE: change this function whenever attributes change
  p | type_;
  p | value_accessed_;
  if (type_ == parameter_integer) {
    p | value_integer_;
  } else if (type_ == parameter_float) {
    p | value_float_;
  } else if (type_ == parameter_logical) {
    p | value_logical_; 
  } else if (type_ == parameter_string) {
    int n = 0;
    if (!up) {
      n=strlen(value_string_);
    }
    p | n;
    // +1 to include '\0' terminator
    if (up) value_string_ = new char [n+1];
    PUParray(p,value_string_,n+1);
  } else if (type_ == parameter_list) {
    int n = 0;
    if (! up) {
      n = value_list_->size();
    }
    p | n;
    if (up) {
      value_list_ = new list_type;
      value_list_->resize(n);
    }
    for (int i=0; i<n; i++) {
      int l_param = 0;
      if (! up) {
	l_param = ((*value_list_)[i] != NULL);
      }
      p | l_param;
      if (l_param) {
	if (up)
	  (*value_list_)[i] = new Param;
	p | *(*value_list_)[i];
      }
    }
  } else if (type_ == parameter_logical_expr) {
    pup_expr_(p,&value_expr_);
  } else if (type_ == parameter_float_expr) {
    pup_expr_(p,&value_expr_);
  } else if (type_ == parameter_unknown) {
    WARNING("Param::pup","parameter type is unknown");
  }
}

//----------------------------------------------------------------------

void Param::pup_expr_ (PUP::er &p, struct node_expr ** node) {
  int l_expr = 0;
  const bool up = p.isUnpacking();
  if (! up ) {
    l_expr = ((*node) != NULL);
  }
  p | l_expr;
  if (!l_expr && up) (*node) = NULL;
  if (l_expr) {

    // PUP type
    if (up) {
      (*node) = new struct node_expr;
    }
    p | (*node)->type;

    // PUP function_name
    int n = 0;
    if (!up && (*node)->function_name) {
      n = strlen((*node)->function_name);
    }
    p | n;
    if (n > 0) {
      if (up) {
	(*node)->function_name = new char [ n + 1];
      }
      PUParray (p,(*node)->function_name,n+1);
    } else {
      if (up) (*node)->function_name = NULL;
    }

    // PUP value (note function_name needs to be before value)
    switch ((*node)->type) {
    case enum_node_operation:
      p | (*node)->op_value;
      break;
    case enum_node_float:
      p | (*node)->float_value;
      break;
    case enum_node_integer:
      p | (*node)->integer_value;
      break;
    case enum_node_variable:
      p | (*node)->var_value;
      break;
    case enum_node_function:
      if (up) {
      	(*node)->fun_value = function_map.at((*node)->function_name);
      }
      break;
    default:
      WARNING("Param::pup_expr_()",
	      "Unknown node type");
      break;
    }

    // Recurse on left and right subexpressions
    pup_expr_(p,&(*node)->left);
    pup_expr_(p,&(*node)->right);
  }
}

//----------------------------------------------------------------------

void Param::set (struct param_struct * node)
/// @param   node  The node from which to copy the type and value
{

  value_accessed_ = false;

  switch (node->type) {
  case enum_parameter_integer:
    set_integer_(node->integer_value);
    break;
  case enum_parameter_float:
    set_float_(node->float_value);
    break;
  case enum_parameter_string:
    set_string_(node->string_value);
    break;
  case enum_parameter_logical:
    set_logical_(node->logical_value);
    break;
  case enum_parameter_list:
    set_list_(node->list_value);
    break;
  case enum_parameter_float_expr:
    set_float_expr_(node->op_value);
    break;
  case enum_parameter_logical_expr:
    set_logical_expr_(node->op_value);
    break;
  case enum_parameter_unknown:
  case enum_parameter_sentinel:
  case enum_parameter_function:
  case enum_parameter_group:
  case enum_parameter_identifier:
    break;
  }
}

//----------------------------------------------------------------------

void Param::dealloc_() 
///
{ 
  switch (type_) {
  case parameter_string: 
    dealloc_string_(); 
    break;
  case parameter_list:   
    dealloc_list_(value_list_); 
    break;
  case parameter_logical_expr:
  case parameter_float_expr:
    dealloc_node_expr_(value_expr_);
    break;
  case parameter_unknown:
  case parameter_integer:
  case parameter_float:
  case parameter_logical:
    break;
  }
} 

//----------------------------------------------------------------------

void Param::write
(
 FILE *      file_pointer,
 std::string full_parameter,
 int         write_type)
/// @param file_pointer    File pointer to which the parameter is written
/// @param full_parameter  Name of this parameter including groups
{

  // Write the parameter assignment

  size_t i_group = full_parameter.rfind(":");
  std::string parameter = (i_group == std::string::npos) ?
    full_parameter : full_parameter.substr(i_group+1,std::string::npos);

  std::string value = value_to_string(write_type);
  if (write_type == param_write_monitor) {
    std::size_t found;
    while ((found=value.find(",")) != std::string::npos ) {
      value.erase(found,1);
    }
    fprintf (file_pointer,"%s = %s\n",
  	     parameter.c_str(),
  	     value.c_str());
  } else {
    fprintf (file_pointer,"%s = %s;\n",
	     parameter.c_str(),
	     value.c_str());
  }
  
}

//----------------------------------------------------------------------

std::string Param::value_to_string (int type)
{
  std::string string_buffer;
  char char_buffer[MAX_BUFFER_LENGTH];

  const std::string list_begin = (type == param_write_libconfig) ? "( " : "[ ";
  const std::string list_end   = (type == param_write_libconfig) ? " )" : " ]";
  const std::string expr_begin = (type == param_write_libconfig) ? "\"" : "";
  const std::string expr_end   = (type == param_write_libconfig) ? "\"" : "";
  
  switch (type_) {
  case parameter_string: 
    string_buffer = std::string("\"") + value_string_ + "\"";
    break;
  case parameter_list:
    string_buffer = list_begin;
    for (size_t i=0; i<value_list_->size(); i++) {
      if ( i > 0 ) string_buffer += ", ";
      string_buffer += (*value_list_)[i]->value_to_string(type);
    }
    string_buffer += list_end;
    break;
  case parameter_logical_expr:
  case parameter_float_expr:
    sprintf_expression(value_expr_,char_buffer);
    string_buffer = expr_begin + char_buffer + expr_end;
    break;
  case parameter_integer:
    sprintf (char_buffer,"%d",value_integer_);
    string_buffer = char_buffer;
    break;
  case parameter_float:
    // '#' format character forces a decimal point, which is required to
    // differentiate an integer from a float type
    sprintf (char_buffer,FLOAT_FORMAT,value_float_);
    string_buffer =  char_buffer;
    break;
  case parameter_logical:
    string_buffer = value_logical_ ? "true" : "false";
    break;
  case parameter_unknown:
    string_buffer = "UNKNOWN";
    break;
  }  
  return string_buffer;
}

//----------------------------------------------------------------------

void Param::dealloc_list_ (list_type * value)
/// @param value List to be deallocated
{
  for (unsigned i=0; i<(*value).size(); i++) {
    if ((*value)[i]->type_ == parameter_list) {
      dealloc_list_ ((*value)[i]->value_list_);
    } else {
      delete ( (*value)[i] );
    }
  }
  delete value;
}

//----------------------------------------------------------------------

void Param::dealloc_node_expr_ (struct node_expr * p)
/// @param p Head node of the tree defining the floating-point 
/// expression to deallocate
{
  if (p->left != NULL)  dealloc_node_expr_(p->left);
  if (p->right != NULL) dealloc_node_expr_(p->right);
  free (p->function_name);
  free (p);
}

//----------------------------------------------------------------------
