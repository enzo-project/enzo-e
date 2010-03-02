/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

/** 
*********************************************************************
*
* @file      param.cpp
* @brief     Helper Param* classes for parameters
* @author    James Bordner
* @date      Sun Oct 11 15:02:08 PDT 2009
* @bug       Probable memory leaks
* @note      
*
* DESCRIPTION 
* 
*    Classes Param* to represent and evaluate various parameter types and
*    expressions
*
* PACKAGES
*
*    NONE
* 
* INCLUDES
*  
*    NONE
*
* PUBLIC FUNCTIONS
*  
* PRIVATE FUCTIONS
*  
* $Id$
*
*********************************************************************
*/

#include "stdio.h"
#include "stdlib.h"

#include "parse.h"
#include "param.hpp"
#include "error.hpp"

      
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

void Param::set (struct param_type * node)
{
  switch (node->type) {
  case enum_parameter_integer:
    set_integer_(node->integer_value);
    break;
  case enum_parameter_scalar:
    set_scalar_(node->scalar_value);
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
  case enum_parameter_scalar_expr:
    set_scalar_expr_(node->op_value);
    break;
  case enum_parameter_logical_expr:
    set_logical_expr_(node->op_value);
    break;
  case enum_parameter_unknown:
  case enum_parameter_sentinel:
  case enum_parameter_identifier:
  case enum_parameter_function:
    break;
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

void Param::dealloc_() { 
  switch (type_) {
  case type_string_: 
    dealloc_string_(); 
    break;
  case type_list_:   
    dealloc_list_(value_list_); 
    break;
  case type_logical_expr_:
  case type_scalar_expr_:
    dealloc_node_expr_(value_expr_);
    break;
  case type_unknown_:
  case type_integer_:
  case type_scalar_:
  case type_logical_:
    break;
  }
} 

// Write parameter to file

void Param::write(FILE * file_pointer,
		  std::string parameter)
{
  // Write the parameter name
  fprintf (file_pointer,"%s ",parameter.c_str());

  // Write the parameter value
  switch (type_) {
  case type_string_: 
    fprintf (file_pointer,"%s\n",value_string_);
    break;
  case type_list_:
    fprintf (file_pointer,"LIST\n");
    INCOMPLETE_MESSAGE("Param::write","Writing lists is not implemented yet");
    break;
  case type_logical_expr_:
  case type_scalar_expr_:
    write_scalar_expr_(file_pointer,value_expr_);
    fprintf (file_pointer,"\n");
    break;
  case type_integer_:
    fprintf (file_pointer,"%d\n",value_integer_);
    break;
  case type_scalar_:
    fprintf (file_pointer,"%g\n",value_scalar_);
    break;
  case type_logical_:
    fprintf (file_pointer,"%s\n",value_logical_ ? "true" : "false");
    break;
  case type_unknown_:
    fprintf (file_pointer,"UNKNOWN\n");
    break;
  }
}


void Param::write_scalar_expr_(FILE * fp,
			       struct node_expr * node)
{
  print_expression(node,fp);
}


// Evaluate scalar expression


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

void Param::evaluate_scalar
(struct node_expr * node, 
 int                n, 
 double *           result, 
 double *           x, 
 double *           y, 
 double *           z, 
 double *           t)
{
  double * left  = NULL;
  double * right = NULL;

  if (node->left) {
    left = new double [n];
    evaluate_scalar(node->left,n,left,x,y,z,t);
  }
  if (node->right) {
    right = new double [n];
    evaluate_scalar(node->right,n,right,x,y,z,t);
  }
      
  int i;
  switch (node->type) {
  case enum_node_operation:
    switch (node->op_value) {
    case enum_op_add: for (i=0; i<n; i++) result[i] = left[i] + right[i]; break;
    case enum_op_sub: for (i=0; i<n; i++) result[i] = left[i] - right[i]; break;
    case enum_op_mul: for (i=0; i<n; i++) result[i] = left[i] * right[i]; break;
    case enum_op_div: for (i=0; i<n; i++) result[i] = left[i] / right[i]; break;
    default:
    case enum_op_le:
    case enum_op_lt:
    case enum_op_ge:
    case enum_op_gt:
    case enum_op_eq:
    case enum_op_ne:
    case enum_op_and:
    case enum_op_or:
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"logical operator %d in scalar expression",
	       node->op_value);
      ERROR_MESSAGE("Param::evaluate_scalar",error_message);
      break;
    }
    break;
  case enum_node_scalar:
    for (i=0; i<n; i++) result[i] = node->scalar_value;
    break;
  case enum_node_integer:
    for (i=0; i<n; i++) result[i] = double(node->integer_value);
    break;
  case enum_node_variable:
    switch (node->var_value) {
    case 'x':	for (i=0; i<n; i++) result[i] = x[i]; break;
    case 'y':	for (i=0; i<n; i++) result[i] = y[i]; break;
    case 'z':	for (i=0; i<n; i++) result[i] = z[i]; break;
    case 't':	for (i=0; i<n; i++) result[i] = t[i]; break;
    default:
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"unknown variable %c in scalar expression",
	       node->var_value);
      ERROR_MESSAGE("Param::evaluate_scalar",error_message);
      break;
    }
    break;
  case enum_node_function:
    for (i=0; i<n; i++) result[i] = (*(node->fun_value))(left[i]);
    break;
  case enum_node_unknown:
  default:
    char error_message[ERROR_MESSAGE_LENGTH];
    sprintf (error_message,"unknown expression type %d",node->type);
    ERROR_MESSAGE("Param::evaluate_scalar",error_message);
    break;
  }

  delete [] left;
  delete [] right;
}

// Evaluate logical expression


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

void Param::evaluate_logical
(struct node_expr * node, 
 int                n, 
 bool   *           result, 
 double *           x, 
 double *           y, 
 double *           z, 
 double *           t)
{
  double * left_scalar  = NULL;
  double * right_scalar = NULL;
  bool * left_logical  = NULL;
  bool * right_logical = NULL;

  // Recurse on left subtree

  if (node->left && node->left->type == enum_node_operation) {
    // left node is an operation
    if (node->op_value == enum_op_and || 
	node->op_value == enum_op_or) {
      // left node is a logical operation
      left_logical = new bool [n];
      evaluate_logical(node->left,n,left_logical,x,y,z,t);
    } else {
      // left node is a scalar operation
      left_scalar = new double [n];
      evaluate_scalar(node->left,n,left_scalar,x,y,z,t);
    }
  } else {
    // left node is a scalar operation
    left_scalar = new double [n];
    evaluate_scalar(node->left,n,left_scalar,x,y,z,t);
  }

  // Recurse on left subtree

  if (node->right && node->right->type == enum_node_operation) {
    // right node exists
    // right node is an operation
    if (node->op_value == enum_op_and || 
	node->op_value == enum_op_or) {
      // right node is a logical operation
      right_logical = new bool [n];
      evaluate_logical(node->right,n,right_logical,x,y,z,t);
    } else {
      // right node is a scalar operation
      right_scalar = new double [n];
      evaluate_scalar(node->right,n,right_scalar,x,y,z,t);
    }
  } else {
    // right node is a scalar operation
    right_scalar = new double [n];
    evaluate_scalar(node->right,n,right_scalar,x,y,z,t);
  }
      
  int i;
  if (node->type == enum_node_operation) {
    switch (node->op_value) {
    case enum_op_le:
      for (i=0; i<n; i++) result[i] = left_scalar[i] <= right_scalar[i];
      break;
    case enum_op_lt:
      for (i=0; i<n; i++) result[i] = left_scalar[i] <  right_scalar[i];
      break;
    case enum_op_ge:
      for (i=0; i<n; i++) result[i] = left_scalar[i] >= right_scalar[i];
      break;
    case enum_op_gt:
      for (i=0; i<n; i++) result[i] = left_scalar[i] >  right_scalar[i];
      break;
    case enum_op_eq:
      // warning: comparing equality of doubles
      for (i=0; i<n; i++) result[i] = left_scalar[i] == right_scalar[i];
      break;
    case enum_op_ne:
      // warning: comparing inequality of doubles
      for (i=0; i<n; i++) result[i] = left_scalar[i] != right_scalar[i];
      break;
    case enum_op_and:
      for (i=0; i<n; i++) result[i] = left_logical[i] && right_logical[i];
      break;
    case enum_op_or:
      for (i=0; i<n; i++) result[i] = left_logical[i] || right_logical[i];
      break;
    default:
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"unknown expression type %d",node->type);
      ERROR_MESSAGE("Param::evaluate_logical",error_message);
      break;
    }
  }

  delete [] left_scalar;
  delete [] right_scalar;
  delete [] left_logical;
  delete [] right_logical;
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

void Param::dealloc_list_ (list_type * value)
{
  for (unsigned i=0; i<(*value).size(); i++) {
    if ((*value)[i]->type_ == type_list_) {
      dealloc_list_ ((*value)[i]->value_list_);
    } else {
      delete ( (*value)[i] );
    }
  }
  delete value;
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

void Param::dealloc_node_expr_ (struct node_expr * p)
{
  if (p->left != NULL)  dealloc_node_expr_(p->left);
  if (p->right != NULL) dealloc_node_expr_(p->right);
  free (p->function_name);
  free (p);
};
