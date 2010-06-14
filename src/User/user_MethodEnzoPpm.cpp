// $Id: user_MethodEnzoPpm.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the MethodEnzoPpm class

#include "user_MethodEnzoPpm.hpp"
#include "error.hpp"
#include "parameters.hpp"

void MethodEnzoPpm::initialize() throw()
{
  // Register method name

  method_name_ = "ppm";

  // Specify arguments

  add_argument_(argument_field, "density",        access_read_write);
  add_argument_(argument_field, "energy_total",   access_read_write);
  add_argument_(argument_field, "energy_internal",access_read_write);
  add_argument_(argument_field, "velocity_x",     access_read_write);
  add_argument_(argument_field, "velocity_y",     access_read_write);
  add_argument_(argument_field, "velocity_z",     access_read_write);

  // Initialize from parameters

  Parameters * parameters = Parameters::instance();

  parameters->set_current_group ("Method","ppm");

  diffusion_  = parameters->value_logical ("diffusion", "false");
  steepening_ = parameters->value_logical ("steepening", "false");
  flattening_ = parameters->value_logical ("flattening", "false");
}

void MethodEnzoPpm::advance_block() throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpm::advance_block","");
}

void MethodEnzoPpm::refresh_face() throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpm::refresh_face","");
}

