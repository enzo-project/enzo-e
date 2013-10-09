// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpml class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpml::EnzoMethodPpml() : Method()
{
}

//----------------------------------------------------------------------

void EnzoMethodPpml::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodPpml::compute_block
( FieldDescr * field_descr,  CommBlock * block ) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveMHDEquations ( field_descr, block->dt() );
}

//----------------------------------------------------------------------
