// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlockFaces.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file field_FieldBlockFaces.cpp
///
/// Detailed description of file field_FieldBlockFaces.cpp

#include "cello.hpp"
#include "field.hpp"

FieldBlockFaces::FieldBlockFaces() throw ()
{
  for (int i=0; i<6; i++) {
    faces_[i] = 0;
  }
}

//----------------------------------------------------------------------

FieldBlockFaces::~FieldBlockFaces() throw ()
{
  for (int i=0; i<6; i++) delete [] faces_[i];
}

//----------------------------------------------------------------------

FieldBlockFaces::FieldBlockFaces(const FieldBlockFaces & field_block_faces) throw ()
/// @param     FieldBlockFaces  Object being copied
{
  INCOMPLETE_MESSAGE("FieldBlockFaces::FieldBlockFaces","");
}

//----------------------------------------------------------------------

FieldBlockFaces & FieldBlockFaces::operator= (const FieldBlockFaces & field_block_faces) throw ()
/// @param     FieldBlockFaces  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE_MESSAGE("FieldBlockFaces::operator =","");
  return *this;
}

//======================================================================

void FieldBlockFaces::copy_from_block()
{
  INCOMPLETE_MESSAGE("FieldBlockFaces::copy_from_block","");
}

//----------------------------------------------------------------------

void FieldBlockFaces::copy_to_block()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::send_begin()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::send_end()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::recv_begin()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::recv_end()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::exchange_begin()
{
}

//----------------------------------------------------------------------

void FieldBlockFaces::exchange_end()
{
}

//----------------------------------------------------------------------
