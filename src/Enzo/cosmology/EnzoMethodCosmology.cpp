// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCosmology.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    Implementation of the EnzoMethodCosmology class

#include "cello.hpp"
#include "charm_simulation.hpp"
#include "enzo.hpp"

// #define DEBUG_COSMO

//----------------------------------------------------------------------

EnzoMethodCosmology::EnzoMethodCosmology() throw()
: Method ()
{
}

//----------------------------------------------------------------------

void EnzoMethodCosmology::compute(Block * block) throw()
{
  auto cosmology = enzo::cosmology();

#ifdef DEBUG_COSMO  
  cosmology->print();
#endif  

  // Monitor current redshift
  Monitor * monitor = cello::monitor();
  if (block->index().is_root()) {
    monitor->print("Method", "%s redshift %.8f",
		   this->name().c_str(),
		   cosmology->current_redshift());
  }
  
  block->compute_done(); 
}
