// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Value.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-01 17:10:11
/// @brief    

#include "problem.hpp"

//======================================================================

Value::Value(Parameters * parameters,
	     const std::string parameter_name) throw()
{
  const int param_type = parameters->type(parameter_name);

  if (param_type == parameter_float_expr || 
      param_type == parameter_float ) {

    Param * param = parameters->param(parameter_name);

    ScalarExpr * scalar_expr = new ScalarExpr(param);
    scalar_expr_list_.push_back(scalar_expr);

    mask_list_.push_back(NULL);

  }  else if (param_type == parameter_list) {

    const int list_length = parameters->list_length(parameter_name);
    for (int index=0; index<list_length; index+=2) {

      Param * param = parameters->param(parameter_name,index);
      ScalarExpr * scalar_expr = new ScalarExpr(param);
     
      scalar_expr_list_.push_back(scalar_expr);

      if (index+1 < list_length) {
	param = parameters->param(parameter_name,index+1);
	Mask * mask = Mask::create(param,parameters);
	mask_list_.push_back(mask);
      } else {
	mask_list_.push_back(NULL);
      }
    }
  } else {
    ERROR2("Value::Value",
	   "Parameter %s is of incorrect type %d",
	   parameter_name.c_str(),param_type);
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

void Value::copy_(const Value & value) throw()
{
  mask_list_.resize(value.mask_list_.size());
  for (size_t i = 0; i<mask_list_.size(); i++) {
    mask_list_[i] = value.mask_list_[i]->clone();
  }
}

//----------------------------------------------------------------------
