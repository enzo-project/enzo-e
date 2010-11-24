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

FieldFaces::FieldFaces(FieldBlock * field_block) throw ()
  : field_block_(field_block)
{

  for (int i=0; i<6; i++) {
    ghost_[i] = 0;
    face_ [i] = 0;
  }
}

//----------------------------------------------------------------------

FieldFaces::~FieldFaces() throw ()
{
  for (int i=0; i<6; i++) delete [] face_[i];
  for (int i=0; i<6; i++) delete [] ghost_[i];
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

void send_init() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::send_init()","");
}
//----------------------------------------------------------------------
void send_begin() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::send_begin()","");
}
//----------------------------------------------------------------------
void send_end() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::send_end()","");
}
//----------------------------------------------------------------------
void send_final() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::send_final()","");
}
	
//----------------------------------------------------------------------
void recv_init() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::recv_init()","");
}
//----------------------------------------------------------------------
void recv_begin() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::recv_begin()","");
}
//----------------------------------------------------------------------
void recv_end() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::recv_end()","");
}
//----------------------------------------------------------------------
void recv_final() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::recv_final()","");
}
	
//----------------------------------------------------------------------
void sendrecv_init() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::sendrecv_init()","");
}
//----------------------------------------------------------------------
void sendrecv_begin() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::sendrecv_begin()","");
}
//----------------------------------------------------------------------
void sendrecv_end() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::sendrecv_end()","");
}
//----------------------------------------------------------------------
void sendrecv_final() throw()
{
  INCOMPLETE_MESSAGE("FieldFaces::sendrecv_final()","");
}

