/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFluxAccretion.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       26 May 2022
/// @brief      Computes accretion rates according to "Flux Accretion" method, as
///             described in Section 5.3 of Bleuler, A & Teyssier, R 2004; MNRAS, 445, 4015-4036

#include "cello.hpp"
#include "enzo.hpp"

//------------------------------------------------------------------

EnzoMethodFluxAccretion::EnzoMethodFluxAccretion
(double accretion_radius_cells,
 double physical_density_threshold_cgs,
 double max_mass_fraction)
  : EnzoMethodAccretion(accretion_radius_cells,
			physical_density_threshold_cgs,
			max_mass_fraction)
{
  // The number of cell widths within the accretion radius cannot be larger than
  // any of the ghost depths minus 1. This is because this method requires computing
  // derivatives, so cells at the edge of the accretion zone must have another cell
  // "outside" it.
  const int * ghost_depth = enzo::config()->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  ASSERT("EnzoMethodFluxAccretion::EnzoMethodFluxAccretion() ",
	 "The accretion radius must be no greater than than the ghost depth minus one."
	 "(ghost depth is 4 cells by default)",
	 accretion_radius_cells_ <= min_ghost_depth - 1);

}

//-------------------------------------------------------------------

void EnzoMethodFluxAccretion::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretion::pup(p); // call parent class pup

  return;
}

//--------------------------------------------------------------------

void EnzoMethodFluxAccretion::compute (Block * block) throw()
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

void EnzoMethodFluxAccretion::compute_(Block * block)
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

    // Also need the field dimensions
    int mx, my, mz;
    field.dimensions (0, &mx, &my, &mz);

    // Loop over batches
    const int nb = particle.num_batches(it);
    for (int ib=0; ib < nb; ib++){

      // Loop over particles in this batch
      const int np = particle.num_particles(it,ib);
      for (int ip=0; ip < np; ip++){

	// Create an EnzoFluxSinkParticle object
	EnzoFluxSinkParticle sp =
	  EnzoFluxSinkParticle(block,ib,ip,accretion_radius,density_threshold);

	// Update properties of sink particle and set values for source fields
	sp.compute(max_mass_fraction_);

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
