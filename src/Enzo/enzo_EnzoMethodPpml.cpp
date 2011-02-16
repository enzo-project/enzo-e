// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpml class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

// void EnzoMethodPpml::initialize(DataDescr * data_descr) throw()
// {
//   // Specify arguments

//   add_argument_(argument_field, "density",    access_read_write);
//   add_argument_(argument_field, "velocity_x", access_read_write);
//   add_argument_(argument_field, "velocity_y", access_read_write);
//   add_argument_(argument_field, "velocity_z", access_read_write);
//   add_argument_(argument_field, "magnetic_x", access_read_write);
//   add_argument_(argument_field, "magnetic_y", access_read_write);
//   add_argument_(argument_field, "magnetic_z", access_read_write);

//   add_argument_(argument_field, "velocity_x_face_x", access_read_write);
//   add_argument_(argument_field, "velocity_y_face_x", access_read_write);
//   add_argument_(argument_field, "velocity_z_face_x", access_read_write);
//   add_argument_(argument_field, "magnetic_x_face_x", access_read_write);
//   add_argument_(argument_field, "magnetic_y_face_x", access_read_write);
//   add_argument_(argument_field, "magnetic_z_face_x", access_read_write);

//   add_argument_(argument_field, "velocity_x_face_y", access_read_write);
//   add_argument_(argument_field, "velocity_y_face_y", access_read_write);
//   add_argument_(argument_field, "velocity_z_face_y", access_read_write);
//   add_argument_(argument_field, "magnetic_x_face_y", access_read_write);
//   add_argument_(argument_field, "magnetic_y_face_y", access_read_write);
//   add_argument_(argument_field, "magnetic_z_face_y", access_read_write);

//   add_argument_(argument_field, "velocity_x_face_z", access_read_write);
//   add_argument_(argument_field, "velocity_y_face_z", access_read_write);
//   add_argument_(argument_field, "velocity_z_face_z", access_read_write);
//   add_argument_(argument_field, "magnetic_x_face_z", access_read_write);
//   add_argument_(argument_field, "magnetic_y_face_z", access_read_write);
//   add_argument_(argument_field, "magnetic_z_face_z", access_read_write);

// }

//----------------------------------------------------------------------

void EnzoMethodPpml::compute_block(DataBlock * data_block,
				   double t,double dt) throw()
{
  INCOMPLETE("EnzoMethodPpml::compute_block","");
}

//----------------------------------------------------------------------
