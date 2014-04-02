// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Value.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "problem.hpp"

//======================================================================

Value::Value(Parameters * parameters,
	     const std::string parameter_name) throw()
{
  if (parameters->type(parameter_name) == parameter_float_expr) {

    Param * param = parameters->param(parameter_name);

    double v,x=1.0,y=2.0,z=3.0,t=7.0;
    param->evaluate_float(1,&v,&x,&y,&z,t);
    ScalarExpr * scalar_expr = new ScalarExpr(param);
    scalar_expr_list_.push_back(scalar_expr);

  }  else if (parameters->type(parameter_name) == parameter_list) {

    const int list_length = parameters->list_length(parameter_name);
    for (int index=0; index<list_length; index+=2) {

      Param * param = parameters->param(parameter_name,index);
      ScalarExpr * scalar_expr = new ScalarExpr(param);
     
      scalar_expr_list_.push_back(scalar_expr);

      if (index+1 < list_length) {
	param = parameters->param(parameter_name,index+1);
	Mask * mask = NULL;
	if ((param->type() == parameter_logical_expr)) {
	  mask = new MaskExpr(param);
	} else if ((param->type() == parameter_string)) {
	  double xm = parameters->list_value_float(0,"Domain:lower",0.0);
	  double xp = parameters->list_value_float(0,"Domain:upper",0.0);
	  double ym = parameters->list_value_float(1,"Domain:lower",0.0);
	  double yp = parameters->list_value_float(1,"Domain:upper",0.0);
	  mask = new MaskPng(param->get_string(),xm,xp,ym,yp);
	} else {
	  ERROR2("Value::Value",
		 "Element %d of list parameter %s is not a logical expression or a string",
		 index+1,parameter_name.c_str());
	}
	mask_list_.push_back(mask);
      } else {
	mask_list_.push_back(NULL);
      }
    }
  } else {
    ERROR1("Value::Value",
	   "Parameter %s is not an expression or a list of expressions",
	   parameter_name.c_str());
  }
}

//----------------------------------------------------------------------

double Value::evaluate (double t, double x, double y, double z) throw ()
{
  double value = 0.0;
  for (int index = (int)scalar_expr_list_.size()-1; index>=0; index--) {
    value = scalar_expr_list_[index]->evaluate(t,x,y,z,mask_list_[index], value);
  }
  return value;
}

//----------------------------------------------------------------------

void Value::evaluate
(double * values, 
 double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z) throw ()
{
  for (int index = (int)scalar_expr_list_.size()-1; index>=0; index--) {
    scalar_expr_list_[index]->evaluate
      (values,t,ndx,nx,x,ndy,ny,y,ndz,nz,z,mask_list_[index], values);
  }
}

//----------------------------------------------------------------------

void Value::copy_(const Value & value) throw()
{
  mask_list_.resize(value.mask_list_.size());
  for (size_t i = 0; i<mask_list_.size(); i++) {
    mask_list_[i] = value.mask_list_[i]->clone();
  }
}

//----------------------------------------------------------------------
