// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPpm class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPpm::EnzoMethodPpm () : Method()
{
  // PPM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoMethodPpm::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

#endif

//----------------------------------------------------------------------

void EnzoMethodPpm::compute_block
(
 FieldDescr * field_descr, CommBlock * block) throw()
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  enzo_block->SolveHydroEquations 
    ( block->cycle(), block->time(), block->dt() );
}

