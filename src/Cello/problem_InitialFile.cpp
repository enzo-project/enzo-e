// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialFile.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @brief    Implementation of the InitialFile class

#include "cello.hpp"

#include "problem.hpp"

//----------------------------------------------------------------------

InitialFile::InitialFile
(int cycle, double time) throw ()
  : Initial (cycle,time)
{
  INCOMPLETE("InitialFile::InitialFile");
}

//----------------------------------------------------------------------

void InitialFile::enforce
(
 const Hierarchy  * hierarchy,
 const FieldDescr * field_descr,
 Block            * block
 ) throw()
{
  INCOMPLETE("InitialFile::enforce");
}

