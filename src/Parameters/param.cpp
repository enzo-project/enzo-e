//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
*********************************************************************
*
* @file      param.cpp
* @brief     Helper Param* classes for parameters
* @author    James Bordner
* @date      Sun Oct 11 15:02:08 PDT 2009
* @bug       
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

void Param::evaluate_scalar
(int                n, 
 double *           result, 
 double *           x, 
 double *           y, 
 double *           z, 
 double *           t)
{
  evaluate_scalar_(value_expr_, n, result, x,y,z,t);
}

// Private evaluate scalar expression


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

void Param::evaluate_scalar_
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
    evaluate_scalar_(node->left,n,left,x,y,z,t);
  }
  if (node->right) {
    right = new double [n];
    evaluate_scalar_(node->right,n,right,x,y,z,t);
  }
      
  int i;
  switch (node->type) {
  case enum_node_operation:
    switch (node->op_value) {
    case enum_op_add:
      for (i=0; i<n; i++) result[i] = left[i] + right[i]; 
      break;
    case enum_op_sub:
      for (i=0; i<n; i++) result[i] = left[i] - right[i]; 
      break;
    case enum_op_mul:
      for (i=0; i<n; i++) result[i] = left[i] * right[i]; 
      break;
    case enum_op_div:
      for (i=0; i<n; i++) result[i] = left[i] / right[i]; 
      break;
    default:
    case enum_op_le:
    case enum_op_lt:
    case enum_op_ge:
    case enum_op_gt:
    case enum_op_eq:
    case enum_op_ne:
    case enum_op_and:
    case enum_op_or:
      sprintf (error_message,"logical operator %d in scalar expression\n",
	       node->op_value);
      ERROR_MESSAGE("Param::evaluate_scalar");
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
      sprintf (error_message,"unknown variable %c in scalar expression\n",
	       node->var_value);
      ERROR_MESSAGE("Param::evaluate_scalar");
      break;
    }
    break;
  case enum_node_function:
    for (i=0; i<n; i++) result[i] = (*(node->fun_value))(left[i]);
    break;
  case enum_node_unknown:
  default:
    sprintf (error_message,"unknown expression type %d\n",node->type);
    ERROR_MESSAGE("Param::evaluate_scalar");
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

// bool Param::evaluate_logical_
// (struct node_expr * node, 
//  int                n, 
//  double *           result, 
//  double *           x, 
//  double *           y, 
//  double *           z, 
//  double *           t)
// {
//   bool result;
//   return result;
// }



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

void Param::dealloc_param_ (struct param_type * p)
{
  struct param_type * head = p;
  struct param_type * q;
  p = p->next;
  while (p->type != enum_parameter_sentinel) {
    q = p->next;
    if (p->type == enum_parameter_list) {
      dealloc_param_ (p);
    } else {
      free ( p );
    }
    p = q;
  }
  free (head);
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
