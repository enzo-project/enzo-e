// See LICENSE_CELLO file for license and copyright information

/// @file     _parameters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:44:21 PDT 2009
/// @brief    Private include file for the \ref Parameters component

#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

/// @enum     parameter_enum
/// @brief    Parameter data type
enum parameter_enum {
  parameter_unknown,
  parameter_integer,
  parameter_float,
  parameter_string,
  parameter_logical,
  parameter_list,
  parameter_float_expr,
  parameter_logical_expr
};

typedef int parameter_type;

#define MAX_BUFFER_LENGTH 10000

extern const char * parameter_type_name [];

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include <limits>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "parse.h"
#include "parameters_Config.hpp"
#include "parameters_Param.hpp"
#include "parameters_ParamNode.hpp"
#include "parameters_Parameters.hpp"

#endif /* _PARAMETERS_HPP */

