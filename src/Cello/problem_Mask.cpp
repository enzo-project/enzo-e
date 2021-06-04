// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Mask.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the Mask class

#include "problem.hpp"

std::shared_ptr<Mask> Mask::create(Parameters * parameters,
				   const std::string &parameter_name,
				   int parameter_index)
{
  Param * param = parameters->param(parameter_name, parameter_index);
  std::shared_ptr<Mask> mask(nullptr);
  if (param) {
    if ((param->type() == parameter_logical_expr)) {
      mask = std::make_shared<MaskExpr> (parameters, parameter_name,
					 parameter_index);
    } else if ((param->type() == parameter_string)) {
      double xm = 0.0;
      double ym = 0.0;
      double xp = 0.0;
      double yp = 0.0;
      if (parameters) {
	xm = parameters->list_value_float(0,"Domain:lower",0.0);
	ym = parameters->list_value_float(1,"Domain:lower",0.0);
	xp = parameters->list_value_float(0,"Domain:upper",0.0);
	yp = parameters->list_value_float(1,"Domain:upper",0.0);
      }
      mask = std::make_shared<MaskPng> (param->get_string(),xm,xp,ym,yp);
    } else {
      ERROR("Mask::create()",
	    "Invalid Mask type: must be logical expression or file name");
    }
  } 
  return mask;

}
//======================================================================

