// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpml class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpml::EnzoMethodPpml( Parameters * parameters )
  : Method(parameters)
{
}

//----------------------------------------------------------------------

void EnzoMethodPpml::compute_block
(
 FieldDescr * field_descr, 
 Block * block,
 int cycle,
 double time,
 double dt
 ) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveMHDEquations ( field_descr, dt);
}

//----------------------------------------------------------------------
