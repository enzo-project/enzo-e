// $Id: user_MethodEnzoPpm.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the MethodEnzoPpm class

#include "data.hpp"
#include "error.hpp"
#include "global.hpp"
#include "user_MethodEnzoPpm.hpp"

#include "cello_hydro.h"

void print_fields(int mx,int my);

//----------------------------------------------------------------------

void MethodEnzoPpm::initialize (DataDescr * data_descr) throw()
{

  // Specify arguments

  add_argument_(argument_field, "density", access_read_write, data_descr);
  add_argument_(argument_field, "total_energy",	access_read_write, data_descr);
  add_argument_(argument_field, "internal_energy",
		access_read_write, data_descr);

  // (get GridRank to only add required velocity fields)

  Parameters * parameters = global_->parameters();
  parameters->set_current_group("Physics");
  GridRank = parameters->value_integer ("dimensions",0);

  if (GridRank >= 1) {
    add_argument_(argument_field, "velocity_x", access_read_write, data_descr);
  }
  if (GridRank >= 2) {  
    add_argument_(argument_field, "velocity_y", access_read_write, data_descr);
  }
  if (GridRank >= 3) {
    add_argument_(argument_field, "velocity_z", access_read_write, data_descr);
  }

  parameters->set_current_group("Method","ppm");
  PPMFlatteningParameter = parameters->value_logical("flattening",true);
  PPMDiffusionParameter  = parameters->value_logical("diffusion",true);
  PPMSteepeningParameter = parameters->value_logical("steepening",true);
}

//----------------------------------------------------------------------

void MethodEnzoPpm::finalize ( DataDescr * data_descr ) throw ()
{
}

//----------------------------------------------------------------------

void MethodEnzoPpm::initialize_block ( DataBlock * data_block ) throw ()
{

}

//----------------------------------------------------------------------

void MethodEnzoPpm::finalize_block ( DataBlock * data_block ) throw ()
{
}

//----------------------------------------------------------------------

void MethodEnzoPpm::advance_block
(
 DataBlock * data_block,
 double t,
 double dt
 ) throw()
{

  SolveHydroEquations (CycleNumber, dt);
}

//----------------------------------------------------------------------

void MethodEnzoPpm::refresh_face() throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpm::refresh_face","");
}

