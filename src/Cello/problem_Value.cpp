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

    ASSERT1("Value::Value()",
	    "param = NULL for parameter %s",
	    parameter_name.c_str(),
	    param != NULL);
	    
    ScalarExpr * scalar_expr = new ScalarExpr(param);
    scalar_expr_list_.push_back(scalar_expr);

    mask_list_.push_back(nullptr);

  }  else if (param_type == parameter_list) {

    const int list_length = parameters->list_length(parameter_name);
    for (int index=0; index<list_length; index+=2) {

      Param * param = parameters->param(parameter_name,index);

      ASSERT2("Value::Value",
	      "param = NULL for parameter %s index %d",
	      parameter_name.c_str(),index,
	      param != NULL);
      
      ScalarExpr * scalar_expr = new ScalarExpr(param);

      scalar_expr_list_.push_back(scalar_expr);

      if (index+1 < list_length) {
	param = parameters->param(parameter_name,index+1);
	mask_list_.push_back(Mask::create(param,parameters));
      } else {
	mask_list_.push_back(nullptr);
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

template <class T>
void Value::evaluate
(T * values, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z) throw ()
{
  for (int index = (int)scalar_expr_list_.size()-1; index>=0; index--) {
    scalar_expr_list_[index]->evaluate
      (values,t,ndx,nx,x,ndy,ny,y,ndz,nz,z,mask_list_[index], values);
  }
}

template void Value::evaluate
(float * values, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z) throw ();

template void Value::evaluate
(double * values, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z) throw ();

template void Value::evaluate
(long double * values, double t,
 int ndx, int nx, double * x,
 int ndy, int ny, double * y,
 int ndz, int nz, double * z) throw ();

//----------------------------------------------------------------------

void Value::copy_(const Value & value) throw()
{
  mask_list_.resize(value.mask_list_.size());
  for (size_t i = 0; i<mask_list_.size(); i++) {
    mask_list_[i] = value.mask_list_[i]->make_clone();
  }
}

//----------------------------------------------------------------------
