// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Expression.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Jun 7 Sun 2021
/// @brief    Implementation of the Expression class

#include "cello.hpp"

#include "parameters.hpp"

//----------------------------------------------------------------------

void Expression::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  const bool up = p.isUnpacking();
  bool is_initialized;
  if (!up){
    is_initialized = initialized();
  }

  p|is_initialized;

  if (is_initialized){
    p|param_name_;
    p|param_index_;
    if (up){
      const Parameters* parameters = cello::parameters();
      (*this) = parameters->construct_or_rebuild_Expression_(param_name_,
							     param_index_);
    }
  }
}

//----------------------------------------------------------------------

void Expression::evaluate_float
( int                      n,
  double *                 result,
  double *                 x,
  double *                 y,
  double *                 z,
  double                   t,
  const struct node_expr * node) const noexcept
/// @param node Head node of the tree defining the floating-point expression
/// @param n Length of the result buffer
/// @param result Array in which to store the expression evaluations
/// @param x Array of X spatial values
/// @param y Array of Y spatial values
/// @param z Array of Z spatial values
/// @param t time value
{
  if (node == nullptr){
    ASSERT("Expression::evaluate_float",
	   "The Expression object has not been initialized.",
	   initialized());
    ASSERT("Expression::evaluate_float",
	   "The Expression object doesn't represent a scalar expression.",
	   is_scalar_expression());
    node = value_expr_;
  }

  double * left  = nullptr;
  double * right = nullptr;

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
           "node is nullptr",
           (node != nullptr));
    ASSERT3("Param::evaluate_float()",
	    "Error in operation %d: left %p right %p",
	    node ? node->op_value:-1,left,right,
	    ((left != nullptr) && (right != nullptr)));
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
	    node?node->fun_value : nullptr, left,
	    left != nullptr);
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

void Expression::evaluate_logical
( int                      n,
  bool *                   result,
  double *                 x,
  double *                 y,
  double *                 z,
  double                   t,
  const struct node_expr * node) const noexcept
/// @param node Head node of the tree defining the floating-point expression
/// @param n Length of the result buffer
/// @param result Array in which to store the expression evaluations
/// @param x Array of X spatial values
/// @param y Array of Y spatial values
/// @param z Array of Z spatial values
/// @param t Array of time values
{
  if (node == nullptr){
    ASSERT("Expression::evaluate_logical",
	   "The Expression object has not been initialized.",
	   initialized());
    ASSERT("Expression::evaluate_logical",
	   "The Expression object doesn't represent a logical expression.",
	   is_logical_expression());
    node = value_expr_;
  }


  double * left_float  = nullptr;
  double * right_float = nullptr;
  bool * left_logical  = nullptr;
  bool * right_logical = nullptr;

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
	    (left_float != nullptr && right_float != nullptr));
    ASSERT3("Param::evaluate_logical()",
	    "Error in operation %d: left_logical %p right_logical %p",
	    op,left_logical,right_logical,
	    (op!=enum_op_and && op!=enum_op_or) ||
	    (left_logical != nullptr && right_logical != nullptr));
	    
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
