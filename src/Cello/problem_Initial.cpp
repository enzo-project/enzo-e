// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Problem class

#include "cello.hpp"
#include "main.hpp"
#include "simulation.hpp"
#include "field.hpp"
#include "problem.hpp"

//----------------------------------------------------------------------

extern CProxy_Simulation  proxy_simulation;

//----------------------------------------------------------------------

void Initial::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  p | cycle_;
  p | time_;

}


//----------------------------------------------------------------------

void Initial::enforce_simulation_
( Simulation * simulation ) throw()
{
  enforce_hierarchy (simulation->hierarchy(), simulation->field_descr());
}

//----------------------------------------------------------------------

void Initial::enforce_hierarchy_
(
 Hierarchy * hierarchy,
 const FieldDescr * field_descr
 ) throw()
{

}

//----------------------------------------------------------------------

void Initial::enforce_block_
(
 Block            * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()
{
  // Enforce fields

  // ItFieldRange it_field (field_descr->field_count());

  // for (it_field.first(); ! it_field.done(); it_field.next()  ) {
  //   const FieldData * field_data = block->field_data();
  //   enforce (field_data,  field_descr, it_field.value());
  // }
  ERROR("Initial::enforce_block_",
	"This function should not be called");
}
