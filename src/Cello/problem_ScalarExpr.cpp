// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-31
/// @brief    Implementation of the ScalarExpr class

#include <cstring>

#include "problem.hpp"

//----------------------------------------------------------------------

ScalarExpr::ScalarExpr
(Param * param) throw()
  : param_(param),
    value_(0)
{
  if (param_->type() == parameter_float) {
    value_ = param_->get_float();
    param_ = 0;
  }
}

//----------------------------------------------------------------------

void ScalarExpr::copy_(const ScalarExpr & scalar_expr) throw()
{
 
  param_ = (scalar_expr.param_) ? new Param(*scalar_expr.param_) : 0;
  value_ = scalar_expr.value_;
}

//----------------------------------------------------------------------

double ScalarExpr::evaluate (double t, double x, double y, double z, 
			     Mask * mask, double deflt) const
{
  double value;
  bool m = mask ? mask->evaluate(t,x,y,z) : true;
  if (m) {
    if (param_)
      param_->evaluate_float(1,&value,&x,&y,&z,t);
    else
      value = value_;
  } else {
    value = deflt;
  }
  return value;
}

//----------------------------------------------------------------------

