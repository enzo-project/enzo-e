// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCosmology.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    Implementation of the EnzoMethodCosmology class

#include "charm_simulation.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodCosmology::EnzoMethodCosmology(const FieldDescr * field_descr) throw()
: Method ()
{
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields();
}

//----------------------------------------------------------------------

void EnzoMethodCosmology::compute(Block * block) throw()
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
  EnzoPhysicsCosmology * cosmology = units->cosmology();
  
  cosmology->print();
  
  block->compute_done(); 
}
