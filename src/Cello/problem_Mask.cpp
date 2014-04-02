// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Mask.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the Mask class

#include "problem.hpp"

Mask * Mask::create(Param * param, Parameters * parameters)
{
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
    ERROR("Mask::create()",
	   "Invalid Mask type: must be logical expression or file name");
  }
  return mask;

}
//======================================================================

