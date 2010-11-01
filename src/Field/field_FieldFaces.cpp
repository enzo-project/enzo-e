// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFaces.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file field_FieldFaces.cpp
///
/// Detailed description of file field_FieldFaces.cpp

#include "cello.hpp"
#include "field.hpp"

FieldFaces::FieldFaces() throw ()
{
  for (int i=0; i<6; i++) {
    faces_[i] = 0;
  }
}

//----------------------------------------------------------------------

FieldFaces::~FieldFaces() throw ()
{
  for (int i=0; i<6; i++) delete [] faces_[i];
}

//----------------------------------------------------------------------

FieldFaces::FieldFaces(const FieldFaces & field_block_faces) throw ()
/// @param     FieldFaces  Object being copied
{
  INCOMPLETE_MESSAGE("FieldFaces::FieldFaces","");
}

//----------------------------------------------------------------------

FieldFaces & FieldFaces::operator= (const FieldFaces & field_block_faces) throw ()
/// @param     FieldFaces  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE_MESSAGE("FieldFaces::operator =","");
  return *this;
}

//======================================================================

void FieldFaces::copy_from_block()
{
  INCOMPLETE_MESSAGE("FieldFaces::copy_from_block","");
}

//----------------------------------------------------------------------

void FieldFaces::copy_to_block()
{
}

//----------------------------------------------------------------------

void FieldFaces::send_begin()
{
}

//----------------------------------------------------------------------

void FieldFaces::send_end()
{
}

//----------------------------------------------------------------------

void FieldFaces::recv_begin()
{
}

//----------------------------------------------------------------------

void FieldFaces::recv_end()
{
}

//----------------------------------------------------------------------

void FieldFaces::exchange_begin()
{
}

//----------------------------------------------------------------------

void FieldFaces::exchange_end()
{
}

//----------------------------------------------------------------------
