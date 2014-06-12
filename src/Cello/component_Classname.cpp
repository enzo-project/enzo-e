// See LICENSE_CELLO file for license and copyright information

/// @file     component_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "component.hpp"

//----------------------------------------------------------------------

Classname::Classname() throw ()
{
  INCOMPLETE("Classname::Classname");
}

//----------------------------------------------------------------------

Classname::Classname(const Classname & classname) throw ()
/// @param     classname  Object being copied
{
  INCOMPLETE("Classname::Classname(Classname)");
}

//----------------------------------------------------------------------

Classname & Classname::operator= (const Classname & classname) throw ()
/// @param     classname  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("Classname::operator=");
  return *this;
}

//----------------------------------------------------------------------

Classname::~Classname() throw ()
{
  INCOMPLETE("Classname::~Classname");
}

/// CHARM++ Pack / Unpack function
void Classname::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

//======================================================================

