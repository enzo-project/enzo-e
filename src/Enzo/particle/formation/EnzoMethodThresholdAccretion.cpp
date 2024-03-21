// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodThresholdAccretion.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       10 March 2022
/// @brief      Implementation of EnzoMethodThresholdAccretion, a class
///             from EnzoMethodAccretion.
///             This method reduces the gas density in the accretion zone around
///             a sink particle to
///             max(density_threshold_,(1-max_mass_fraction)*density),
///             and adds mass and momentum lost by the gas to the sink particle.

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

//-------------------------------------------------------------------------------------------

EnzoMethodThresholdAccretion::EnzoMethodThresholdAccretion
(double accretion_radius_cells,
 double physical_density_threshold_cgs,
 double max_mass_fraction)
  : EnzoMethodAccretion(accretion_radius_cells,
			physical_density_threshold_cgs,
			max_mass_fraction)
{

}

//---------------------------------------------------------------------------------------------

void EnzoMethodThresholdAccretion::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretion::pup(p); // call parent class pup

  return;
}

//---------------------------------------------------------------------------------------------

void EnzoMethodThresholdAccretion::compute (Block * block) throw()
{

  if (cello::is_initial_cycle(InitCycleKind::fresh_or_noncharm_restart)) {
    do_checks_(block);
  }

  if (block->is_leaf()) {
    this->compute_(block);
  }
  else block->compute_done();

  return;
}

//----------------------------------------------------------------------------------------------

void EnzoMethodThresholdAccretion::compute_(Block * block)

{
  // Only need to do anything if there are sink particles on this block.
  Particle particle = block->data()->particle();
  int it = particle.type_index("sink");
  int num_particles = particle.num_particles(it);

  // Get density threshold in code units for this cycle (value will change in
  // cosmological simultions.
  const double density_threshold =
    physical_density_threshold_cgs_ / enzo::units()->density();

  if (num_particles > 0) {

    // Get pointer to density field data
    Field field = block->data()->field();
    enzo_float * density = (enzo_float*) field.values("density");

    // Set accretion radius to be accretion_radius_cells_ multipled by minimum cell width
    double hx, hy, hz;
    block->cell_width(&hx, &hy, &hz);
    const double min_cell_width = std::min(hx,std::min(hy,hz));
    const double accretion_radius = accretion_radius_cells_ * min_cell_width;

    // Also need the field dimensions
    int mx, my, mz;
    field.dimensions (0, &mx, &my, &mz);

    // Loop over batches
    const int nb = particle.num_batches(it);
    for (int ib=0; ib < nb; ib++){

      // Loop over particles in this batch
      const int np = particle.num_particles(it,ib);
      for (int ip=0; ip < np; ip++){

	// Create an EnzoSinkParticle object
	EnzoSinkParticle sp =  EnzoSinkParticle(block,ib,ip,accretion_radius);

	// Loop over all cells which contain the accretion zone
	for (int k = sp.min_ind_z(); k <= sp.max_ind_z(); k++){
	  for (int j = sp.min_ind_y(); j <= sp.max_ind_y(); j++){
	    for (int i = sp.min_ind_x(); i <= sp.max_ind_x(); i++){

	      // Check if cell is in accretion zone
	      if (sp.cell_in_accretion_zone(i, j, k)){
		const int index = INDEX(i,j,k,mx,my);
		if (density[index] > density_threshold){

		  const enzo_float density_change =
		    std::min(density[index] - density_threshold,
			     max_mass_fraction_ * density[index]);

		  // Update sink particle data and source fields due to accretion
		  // from this cell
		  sp.update(density_change, index);

		} // if density is above threshold
	      } // if cell is in accretion zone
	    }
	  }
	} // Triple loop over cells containing accretion zone

	// Write the sink particle data to the particle attribute array
	sp.write_particle_data();

      } // Loop over particles in this batch
    } // Loop over batches
  } // if (num_particles > 0)

  // Start the refresh
  EnzoBlock * enzo_block = enzo::block(block);
  cello::refresh(ir_accretion_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_start(ir_accretion_, CkIndex_EnzoBlock::p_method_accretion_end());

  return;
}
