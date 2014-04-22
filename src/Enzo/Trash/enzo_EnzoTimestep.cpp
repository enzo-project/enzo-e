// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestep.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoTimestep class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoTimestep::EnzoTimestep () throw()
  : Timestep()
{

}

//----------------------------------------------------------------------

void EnzoTimestep::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
  Timestep::pup(p);
}

//----------------------------------------------------------------------

double EnzoTimestep::evaluate ( const FieldDescr * field_descr,
			       CommBlock * comm_block ) throw()
{

}

