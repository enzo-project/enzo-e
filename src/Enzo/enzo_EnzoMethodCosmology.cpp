// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCosmology.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    Implementation of the EnzoMethodCosmology class

#include "charm_simulation.hpp"
#include "enzo.hpp"

// #define DEBUG_COSMO

//----------------------------------------------------------------------

EnzoMethodCosmology::EnzoMethodCosmology(const FieldDescr * field_descr) throw()
: Method ()
{
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     enzo_sync_id_method_cosmology);
  //  refresh(ir)->add_all_fields();
}

//----------------------------------------------------------------------

void EnzoMethodCosmology::compute(Block * block) throw()
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
  EnzoPhysicsCosmology * cosmology = units->cosmology();

#ifdef DEBUG_COSMO  
  cosmology->print();
#endif  

  // Monitor current redshift
  Monitor * monitor = block->simulation()->monitor();
  if (block->index().is_root()) {
    monitor->print("Method", "%s redshift %.8f",
		   this->name().c_str(),
		   cosmology->current_redshift());
  }
  
  block->compute_done(); 
}
