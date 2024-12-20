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
    sprintf_expression(value_expr_,char_buffer,sizeof(char_buffer));
    string_buffer = expr_begin + char_buffer + expr_end;
    break;
  case parameter_integer:
    snprintf (char_buffer,sizeof(char_buffer),"%d",value_integer_);
    string_buffer = char_buffer;
    break;
  case parameter_float:
    // '#' format character forces a decimal point, which is required to
    // differentiate an integer from a float type
    snprintf (char_buffer,sizeof(char_buffer),FLOAT_FORMAT,value_float_);
    string_buffer =  char_buffer;
    break;
  case parameter_logical:
    string_buffer = value_logical_ ? "true" : "false";
    break;
  case parameter_unknown:
    string_buffer = std::string("\"") + "UNKNOWN" + "\"";
    break;
  default:
    // if type_ is something different, just set the string to be "dummy".
    string_buffer = std::string("\"") + "dummy" + "\"";
    break;
  }  
  return string_buffer;
}
 
//----------------------------------------------------------------------

void Param::evaluate_float
(int                n, 
 double *           result, 
 double *           x, 
 double *           y, 
 double *           z, 
 double             t,
 struct node_expr * node)
/// @param node Head node of the tree defining the floating-point expression
/// @param n Length of the result buffer
/// @param result Array in which to store the expression evaluations
/// @param x Array of X spatial values
/// @param y Array of Y spatial values
/// @param z Array of Z spatial values
/// @param t time value
{
  if (node == 0) node = value_expr_;

  double * left  = NULL;
  double * right = NULL;

  value_accessed_ = true;

  if (node->left) {
    left = new double [n];
    evaluate_float(n,left,x,y,z,t,node->left);
  }
  if (node->right) {
    right = new double [n];
    evaluate_float(n,right,x,y,z,t,node->right);
  }
  
  int i;
  switch (node->type) {
  case enum_node_operation:
    ASSERT("Param::evaluate_float()",
           "node is NULL",
           (node != NULL));
    ASSERT3("Param::evaluate_float()",
	    "Error in operation %d: left %p right %p",
	    node ? node->op_value:-1,left,right,
	    ((left != NULL) && (right != NULL)));
    switch (node->op_value) {
    case enum_op_add: for (i=0; i<n; i++) result[i] = left[i] + right[i]; break;
    case enum_op_sub: for (i=0; i<n; i++) result[i] = left[i] - right[i]; break;
    case enum_op_mul: for (i=0; i<n; i++) result[i] = left[i] * right[i]; break;
    case enum_op_div: for (i=0; i<n; i++) result[i] = left[i] / right[i]; break;
    case enum_op_pow: for (i=0; i<n; i++) result[i] = pow(left[i], right[i]); break;
    default:
    case enum_op_le:
    case enum_op_lt:
    case enum_op_ge:
    case enum_op_gt:
    case enum_op_eq:
    case enum_op_ne:
    case enum_op_and:
    case enum_op_or:
      ERROR1("Param::evaluate_float",
	     "logical operator %d in floating-point expression",
	     node->op_value);
      break;
    }
    break;
  case enum_node_float:
    for (i=0; i<n; i++) result[i] = node->float_value;
    break;
  case enum_node_integer:
    for (i=0; i<n; i++) result[i] = double(node->integer_value);
    break;
  case enum_node_variable:
    switch (node->var_value) {
    case 'x': if (x) for (i=0; i<n; i++) result[i] = x[i]; break;
    case 'y': if (y) for (i=0; i<n; i++) result[i] = y[i]; break;
    case 'z': if (z) for (i=0; i<n; i++) result[i] = z[i]; break;
    case 't': for (i=0; i<n; i++) result[i] = t;    break;
    default:
      ERROR1("Param::evaluate_float",
	     "unknown variable %c in floating-point expression",
	     node->var_value);
      break;
    }
    break;
  case enum_node_function:
    ASSERT2("Param::evaluate_float()",
	    "Error in function %p: left %p",
	    node?node->fun_value : NULL, left,
	    left != NULL);
    for (i=0; i<n; i++) result[i] = (*(node->fun_value))(left[i]);
    break;
  case enum_node_unknown:
  default:
    ERROR1("Param::evaluate_float",
	   "unknown expression type %d",
	   node->type);
    break;
  }

  delete [] left;
  delete [] right;
}

//----------------------------------------------------------------------

void Param::evaluate_logical
(int                n, 
 bool   *           result, 
 double *           x, 
 double *           y, 
 double *           z, 
 double             t,
 struct node_expr * node)
/// @param node Head node of the tree defining the floating-point expression
/// @param n Length of the result buffer
/// @param result Array in which to store the expression evaluations
/// @param x Array of X spatial values
/// @param y Array of Y spatial values
/// @param z Array of Z spatial values
/// @param t Array of time values
{
  if (node == 0) node = value_expr_;

  double * left_float  = NULL;
  double * right_float = NULL;
  bool * left_logical  = NULL;
  bool * right_logical = NULL;

  value_accessed_ = true;

  // Recurse on left subtree

  if (node->left && (node->left->type == enum_node_operation)) {
    // left node is an operation
    if ((node->op_value == enum_op_and) || 
	(node->op_value == enum_op_or)) {
      // left node is a logical operation
      left_logical = new bool [n];
      evaluate_logical(n,left_logical,x,y,z,t,node->left);
    } else {
      // left node is a floating-point operation
      left_float = new double [n];
      evaluate_float(n,left_float,x,y,z,t,node->left);
    }
  } else {
    // left node is a floating-point operation
    left_float = new double [n];
    evaluate_float(n,left_float,x,y,z,t,node->left);
  }

  // Recurse on left subtree

  if (node->right && (node->right->type == enum_node_operation)) {
    // right node exists
    // right node is an operation
    if ((node->op_value == enum_op_and) || 
	(node->op_value == enum_op_or)) {
      // right node is a logical operation
      right_logical = new bool [n];
      evaluate_logical(n,right_logical,x,y,z,t,node->right);
    } else {
      // right node is a floating-point operation
      right_float = new double [n];
      evaluate_float(n,right_float,x,y,z,t,node->right);
    }
  } else {
    // right node is a floating-point operation
    right_float = new double [n];
    evaluate_float(n,right_float,x,y,z,t,node->right);
  }
      
  int i;
  if (node->type == enum_node_operation) {
    int op = node->op_value;
    ASSERT3("Param::evaluate_logical()",
	    "Error in operation %d: left_float %p right_float %p",
	    op,left_float,right_float,
	    (op==enum_op_and || op==enum_op_or) ||
	    (left_float != NULL && right_float != NULL));
    ASSERT3("Param::evaluate_logical()",
	    "Error in operation %d: left_logical %p right_logical %p",
	    op,left_logical,right_logical,
	    (op!=enum_op_and && op!=enum_op_or) ||
	    (left_logical != NULL && right_logical != NULL));
	    
    switch (op) {
    case enum_op_le:
      for (i=0; i<n; i++) result[i] = left_float[i] <= right_float[i];
      break;
    case enum_op_lt:
      for (i=0; i<n; i++) result[i] = left_float[i] <  right_float[i];
      break;
    case enum_op_ge:
      for (i=0; i<n; i++) result[i] = left_float[i] >= right_float[i];
      break;
    case enum_op_gt:
      for (i=0; i<n; i++) result[i] = left_float[i] >  right_float[i];
      break;
    case enum_op_eq:
      // warning: comparing equality of doubles
      for (i=0; i<n; i++) result[i] = left_float[i] == right_float[i];
      break;
    case enum_op_ne:
      // warning: comparing inequality of doubles
      for (i=0; i<n; i++) result[i] = left_float[i] != right_float[i];
      break;
    case enum_op_and:
      for (i=0; i<n; i++) result[i] = left_logical[i] && right_logical[i];
      break;
    case enum_op_or:
      for (i=0; i<n; i++) result[i] = left_logical[i] || right_logical[i];
      break;
    default:
      ERROR1("Param::evaluate_logical",
	     "unknown expression type %d",
	     node->type);
      break;
    }
  }

  delete [] left_float;
  delete [] right_float;
  delete [] left_logical;
  delete [] right_logical;
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
