/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodBondiHoyleAccretion.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @author     John Regan (john.regan@mu.ie)
/// @date       13 April 2022
/// @brief      Computes accretion rates according to Bondi-Hoyle model.
///             See Krumholz+ 2004, ApJ, 611, 399 for details.
///

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

//------------------------------------------------------------------

EnzoMethodBondiHoyleAccretion::EnzoMethodBondiHoyleAccretion
(double accretion_radius_cells,
 double physical_density_threshold_cgs,
 double max_mass_fraction)
  : EnzoMethodAccretion(accretion_radius_cells,
			physical_density_threshold_cgs,
			max_mass_fraction)
{

}

//-------------------------------------------------------------------

void EnzoMethodBondiHoyleAccretion::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretion::pup(p); // call parent class pup

  return;
}

//--------------------------------------------------------------------

void EnzoMethodBondiHoyleAccretion::compute (Block * block) throw()
{

  const auto cycle_simulation = enzo::simulation()->state()->cycle();
  const auto cycle_initial = enzo::config()->initial_cycle;
  if ( cycle_simulation == cycle_initial )
    do_checks_(block);

  if (block->is_leaf()) {
    this->compute_(block);
  }
  else block->compute_done();

  return;
}

//-------------------------------------------------------------------------

void EnzoMethodBondiHoyleAccretion::compute_(Block * block)
{
  // Only need to do anything if there are sink particles on this block.
  Particle particle = block->data()->particle();
  int it = particle.type_index("sink");
  int num_particles = particle.num_particles(it);

  if (num_particles > 0) {

    // Get density threshold in code units for this cycle (value will change in
    // cosmological simultions.
    const double density_threshold =
      physical_density_threshold_cgs_ / enzo::units()->density();

    // Get pointer to density field data
    Field field = block->data()->field();
    enzo_float * density = (enzo_float*) field.values("density");

    // Set accretion radius to be accretion_radius_cells_ multipled by minimum cell width
    double hx, hy, hz;
    block->cell_width(&hx, &hy, &hz);
    const double min_cell_width = std::min(hx,std::min(hy,hz));
    const double accretion_radius = accretion_radius_cells_ * min_cell_width;

    // Get gravitational constant in code units
    const double const_G =
      enzo::grav_constant_cgs() * enzo::units()->density() *
      enzo::units()->time() * enzo::units()->time();

    // Also need the field dimensions
    int mx, my, mz;
    field.dimensions (0, &mx, &my, &mz);

    // Loop over batches
    const int nb = particle.num_batches(it);
    for (int ib=0; ib < nb; ib++){

      // Loop over particles in this batch
      const int np = particle.num_particles(it,ib);
      for (int ip=0; ip < np; ip++){

	// Create an EnzoBondiHoyleSinkParticle object
	EnzoBondiHoyleSinkParticle sp =
	  EnzoBondiHoyleSinkParticle(block,ib,ip,accretion_radius,const_G);

	// Update properties of sink particle and set values for source fields
	sp.compute(density_threshold, max_mass_fraction_);

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
