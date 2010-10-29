// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file     enzo_MethodEnzoPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the MethodEnzoPpml class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "error.hpp"
#include "user.hpp"
#include "enzo.hpp"
#include "parameters.hpp"

//----------------------------------------------------------------------

void MethodEnzoPpml::initialize(DataDescr * data_descr) throw()
{
  // Specify arguments

  add_argument_(argument_field, "density",    access_read_write);
  add_argument_(argument_field, "velocity_x", access_read_write);
  add_argument_(argument_field, "velocity_y", access_read_write);
  add_argument_(argument_field, "velocity_z", access_read_write);
  add_argument_(argument_field, "magnetic_x", access_read_write);
  add_argument_(argument_field, "magnetic_y", access_read_write);
  add_argument_(argument_field, "magnetic_z", access_read_write);

  add_argument_(argument_field, "velocity_x_face_x", access_read_write);
  add_argument_(argument_field, "velocity_y_face_x", access_read_write);
  add_argument_(argument_field, "velocity_z_face_x", access_read_write);
  add_argument_(argument_field, "magnetic_x_face_x", access_read_write);
  add_argument_(argument_field, "magnetic_y_face_x", access_read_write);
  add_argument_(argument_field, "magnetic_z_face_x", access_read_write);

  add_argument_(argument_field, "velocity_x_face_y", access_read_write);
  add_argument_(argument_field, "velocity_y_face_y", access_read_write);
  add_argument_(argument_field, "velocity_z_face_y", access_read_write);
  add_argument_(argument_field, "magnetic_x_face_y", access_read_write);
  add_argument_(argument_field, "magnetic_y_face_y", access_read_write);
  add_argument_(argument_field, "magnetic_z_face_y", access_read_write);

  add_argument_(argument_field, "velocity_x_face_z", access_read_write);
  add_argument_(argument_field, "velocity_y_face_z", access_read_write);
  add_argument_(argument_field, "velocity_z_face_z", access_read_write);
  add_argument_(argument_field, "magnetic_x_face_z", access_read_write);
  add_argument_(argument_field, "magnetic_y_face_z", access_read_write);
  add_argument_(argument_field, "magnetic_z_face_z", access_read_write);

}

//----------------------------------------------------------------------

void MethodEnzoPpml::finalize(DataDescr * data_descr) throw()
{}

//----------------------------------------------------------------------

void MethodEnzoPpml::initialize_block (DataBlock * data_block) throw()
{
}

//----------------------------------------------------------------------

void MethodEnzoPpml::finalize_block (DataBlock * data_block) throw()
{}

//----------------------------------------------------------------------

void MethodEnzoPpml::advance_block(DataBlock * data_block,
				   double t,double dt) throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpml::advance_block","");
}

//----------------------------------------------------------------------
