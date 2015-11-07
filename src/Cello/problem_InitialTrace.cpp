// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-11-06
/// @brief    Implementation of the InitialTrace class

#include "problem.hpp"

//----------------------------------------------------------------------

void InitialTrace::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

//======================================================================

void InitialTrace::enforce_block
 ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw()

{
  Field field (block->data()->field());
  Particle particle (block->data()->particle());
}
