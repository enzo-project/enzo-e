// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Sep 21 16:12:22 PDT 2010
/// @brief    Brief description of file mesh_Mesh.cpp

#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh::Mesh(DataDescr * data_descr) throw ()
  : tree_(NULL)
{

}

//----------------------------------------------------------------------

Mesh::~Mesh() throw ()
{
}

//----------------------------------------------------------------------

Mesh::Mesh(const Mesh & mesh) throw ()
/// @param     mesh  Object being copied
{
}

//----------------------------------------------------------------------

Mesh & Mesh::operator= (const Mesh & mesh) throw ()
/// @param     mesh  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}
