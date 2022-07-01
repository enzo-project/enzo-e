/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretion.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   24 February 2022
/// @brief  Implementation of EnzoMethodAccretion, a base class
///         for "accretion" methods. These methods compute
///         the accretion rate onto sink particles, remove mass,
///         momentum, and energy from the gas, and add mass and
///         momentum to the sink particle.

#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodAccretion::EnzoMethodAccretion
(double accretion_radius_cells,
 double physical_density_threshold_cgs,
 double max_mass_fraction)
  : Method(),
    accretion_radius_cells_(accretion_radius_cells),
    physical_density_threshold_cgs_(physical_density_threshold_cgs),
    max_mass_fraction_(max_mass_fraction),
    ir_accretion_(-1)
{
  // This method requires three dimensions.
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion()",
	 "EnzoMethodAccretion requires that we run a 3D problem (Domain: rank = 3)",
	 cello::rank());

  // Check if Adapt:min_face_rank is 0.
  ASSERT("EnzoMethodAccretion",
	 "Adapt:min_face_rank parameter must be equal to 0.",
	 enzo::config()->adapt_min_face_rank == 0);

  // Check if we are running in unigrid mode (will get rid of this in future)
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion()",
	 "EnzoMethodAccretion requires unigrid mode (Adapt : max_level = 0). "
	 "In future, we may put in a refinement condition that blocks containing "
	 "sink particles are at the highest refinement level.",
	 enzo::config()->mesh_max_level == 0);

  // Check that max_mass_fraction_ is between 0 and 1
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion()",
	 "Method:accretion:max_mass_fraction must be between 0 and 1 inclusive.",
	  max_mass_fraction_ >= 0.0 && max_mass_fraction_ <= 1.0);

  const int * ghost_depth = enzo::config()->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  // The number of cell widths within the accretion radius cannot be larger than
  // any of the ghost depths
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion() ",
	 "The accretion radius must be no greater than the ghost depth"
	 "(4 cells by default)",
	 accretion_radius_cells_ <= min_ghost_depth);

  // Define required fields
  cello::define_field("density_source");
  cello::define_field("density_source_accumulate");
  cello::define_field("mom_dens_x_source");
  cello::define_field("mom_dens_x_source_accumulate");
  cello::define_field("mom_dens_y_source");
  cello::define_field("mom_dens_y_source_accumulate");
  cello::define_field("mom_dens_z_source");
  cello::define_field("mom_dens_z_source_accumulate");

  // Initial refresh: refresh all fields and sink particles
  int it = cello::particle_descr()->type_index("sink");
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
  refresh->add_particle(it);

  // Second refresh: add source fields values in ghost zones to source_accumulate fields
  // values in active cells in neighbouring blocks
  ir_accretion_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_accretion_,name()+":add");
  Refresh * refresh_accretion = cello::refresh(ir_accretion_);

  refresh_accretion->set_accumulate(true);
  refresh_accretion->add_field_src_dst
     ("density_source","density_source_accumulate");
  refresh_accretion->add_field_src_dst
    ("mom_dens_x_source","mom_dens_x_source_accumulate");
  refresh_accretion->add_field_src_dst
    ("mom_dens_y_source","mom_dens_y_source_accumulate");
  refresh_accretion->add_field_src_dst
    ("mom_dens_z_source","mom_dens_z_source_accumulate");
  refresh_accretion->add_field_src_dst
    ("te_dens_source","te_dens_source_accumulate");

  refresh_accretion->set_callback(CkIndex_EnzoBlock::p_method_accretion_end());

}

void EnzoMethodAccretion::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | accretion_radius_cells_;
  p | physical_density_threshold_cgs_;
  p | max_mass_fraction_;
  p | ir_accretion_;

  return;
}


//------------------------------------------------------------------
//   This does nothing - business is done in derived
//   classes
void EnzoMethodAccretion::compute(Block *block) throw() {
  EnzoBlock * enzo_block = enzo::block(block);
  cello::refresh(ir_accretion_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_start(ir_accretion_, CkIndex_EnzoBlock::p_method_accretion_end());
  return;
}

//---------------------------------------------------------------------------------------------
void EnzoBlock::p_method_accretion_end()
{
  EnzoMethodAccretion * method = static_cast<EnzoMethodAccretion*> (this->method());
  method->update_fields(this);
  compute_done();
  return;
}

//--------------------------------------------------------------------------------------------
void EnzoMethodAccretion::update_fields(EnzoBlock * enzo_block) throw()
{

  Field field = enzo_block->data()->field();
  int gx,gy,gz;
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  enzo_float * density = (enzo_float*) field.values("density");
  enzo_float * vx = (enzo_float*) field.values("velocity_x");
  enzo_float * vy = (enzo_float*) field.values("velocity_y");
  enzo_float * vz = (enzo_float*) field.values("velocity_z");
  enzo_float * specific_te = (enzo_float*) field.values("total_energy");

  enzo_float * density_source = (enzo_float*) field.values("density_source");
  enzo_float * density_source_accumulate =
    (enzo_float*) field.values("density_source_accumulate");
  enzo_float * mom_dens_x_source = (enzo_float*) field.values("mom_dens_x_source");
  enzo_float * mom_dens_x_source_accumulate =
    (enzo_float*) field.values("mom_dens_x_source_accumulate");
  enzo_float * mom_dens_y_source = (enzo_float*) field.values("mom_dens_y_source");
  enzo_float * mom_dens_y_source_accumulate =
    (enzo_float*) field.values("mom_dens_y_source_accumulate");
  enzo_float * mom_dens_z_source = (enzo_float*) field.values("mom_dens_z_source");
  enzo_float * mom_dens_z_source_accumulate =
    (enzo_float*) field.values("mom_dens_z_source_accumulate");

  // Loop over all active cells, and update the density, velocity, and
  // energy fields, and magnetic fields if necessary
  for (int iz = gz; iz < mz - gz; iz++){
    for (int iy = gy; iy < my - gy; iy++){
      for (int ix = gx; ix < mx - gx; ix++){

	const int i = INDEX(ix,iy,iz,mx,my);

	// Update density
	const enzo_float old_dens = density[i];
	density[i] += density_source[i] + density_source_accumulate[i];

	// Update color fields
	EnzoMethodStarMaker::rescale_densities(enzo_block,i,density[i] / old_dens);

	// Update velocity
	const enzo_float old_vx = vx[i];
	const enzo_float new_mom_dens_x =
	  old_dens * old_vx + mom_dens_x_source[i] + mom_dens_x_source_accumulate[i];
	vx[i] = new_mom_dens_x / density[i];
	const enzo_float old_vy = vy[i];
	const enzo_float new_mom_dens_y =
	  old_dens * old_vy + mom_dens_y_source[i] + mom_dens_y_source_accumulate[i];
	vy[i] = new_mom_dens_y / density[i];
	const enzo_float old_vz = vz[i];
	const enzo_float new_mom_dens_z =
	  old_dens * old_vz + mom_dens_z_source[i] + mom_dens_z_source_accumulate[i];
	vz[i] = new_mom_dens_z / density[i];

	// Change in specific total energy is just due to change in specific kinetic
	// energy
	const enzo_float old_specific_ke =
	  0.5 * ( old_vx * old_vx + old_vy * old_vy + old_vz * old_vz );
	const enzo_float new_specific_ke =
	  0.5 * ( vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] );
	const enzo_float old_specific_te = specific_te[i];
	specific_te[i] += new_specific_ke - old_specific_ke;

      }
    }
  } // Loop over active cells

  // Loop over all cells (including ghost zones), setting all the "source" and
  // "source accumulate" fields to zero, ready for the next cycle.
  for (int i = 0; i < mx*my*mz; i++){
    density_source[i] = 0.0;
    density_source_accumulate[i] = 0.0;
    mom_dens_x_source[i] = 0.0;
    mom_dens_x_source_accumulate[i] = 0.0;
    mom_dens_y_source[i] = 0.0;
    mom_dens_y_source_accumulate[i] = 0.0;
    mom_dens_z_source[i] = 0.0;
    mom_dens_z_source_accumulate[i] = 0.0;
  }
  return;
}

// Required
double EnzoMethodAccretion::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodAccretion::do_checks_(const Block *block) throw()
{
    // Check if merge_sinks method precedes accretion method
    ASSERT("EnzoMethodAccretion",
	   "merge_sinks must precede accretion_compute",
	   enzo::problem()->method_precedes("merge_sinks",
					    "accretion"));

    // Check if merging radius is at least twice that of the accretion
    // radius
    ASSERT("EnzoMethodAccretion",
	   "Merging radius (Method:merge_sinks:merging_radius_cells "
	   "must be at least twice the accretion radius "
	   "(Method:accretion_compute:accretion_radius).",
	   enzo::config()->method_merge_sinks_merging_radius_cells >=
	   2.0 * accretion_radius_cells_);

    // Check if either PPM or VL+CT method is being used.
    ASSERT("EnzoMethodAccretion",
	   "accretion requires ppm or vlct methods.",
	   enzo::problem()->method_exists("mhd_vlct") ||
	   enzo::problem()->method_exists("ppm"));

    // Check if density threshold is at least as large as the density floor
    //
    // The use of the density_dbl_prec() method is a short term hack to always
    // return the density floor in double precision.
    // TODO: refactor this method to use the density() method instead (which
    // returns the floor in a precision of enzo_float)
    const double density_floor =
      enzo::fluid_props()->fluid_floor_config().density_dbl_prec();

    // Get density threshold in code units.
    // In a cosmological simulation, the density unit is the mean matter density
    // of the universe, which decreases with time, which means that the value of
    // a fixed physical density quantity will increase with time. So if the density
    // threshold is above the density floor at the start of the simulation, it is
    // guaranteed to be abolve the density floor at all times.
    const double density_threshold =
      physical_density_threshold_cgs_ / enzo::units()->density();

    CkPrintf("physical_density_threshold_cgs_ = %.15e\n"
             "density_units = %.15e\n"
             "density_floor = %.15e\n",
             physical_density_threshold_cgs_,
             enzo::units()->density(),
             density_floor);

    ASSERT2("EnzoMethodAccretion",
	    "Density threshold must be at least as large as the density "
	    "floor. At start of simulation, density threshold is %g and "
            "density floor is %g (code density units).",
	    density_threshold,
	    density_floor,
	    density_threshold >= density_floor);

    // The accretion radius must be at least as large as half the
    // diagonal width of a cell, to ensure that at least one cell
    // center is always within the accretion zone.
    // To check this, we compute the diagonal cell width divided by
    // the minimum cell width. accretion_radius_cells_ must be at least
    // half of this value.
    // In addition, since accretion_radius_cells_ must be less that the
    // minimum ghost depth, this value must be less that twice the
    // minimum ghost depth.

    // Get the cell widths
    double hx, hy, hz;
    block->cell_width(&hx, &hy, &hz);
    const double min_cell_width = std::min(hx,std::min(hy,hz));

    // Compute diagonal cell width divided by the minimum cell width,
    // note that this ratio is the same at all refinement levels.
    const double diagonal_over_minimum =
      sqrt(hx * hx + hy * hy + hz * hz) / min_cell_width;

    // Get the minimum ghost depth
    const int * ghost_depth = enzo::config()->field_ghost_depth;
    const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

    // Check that diagonal_over_minimum is less than twice the
    // minimum ghost depth.
    ASSERT2("EnzoMethodAccretion",
	    "The diagonal cell width divided by the minimum "
	    "cell width is %g, and the minimum ghost depth "
	    "is %d. The former must be less than twice the "
	    "latter. It is advised to adjust this by changing "
	    "the dimensions of the root level mesh, via the "
	    "Mesh:root_size parameter.",
	    diagonal_over_minimum, min_ghost_depth,
	    diagonal_over_minimum < 2 * min_ghost_depth);

    // Check that diagonal_over_minimum is less than twice
    // accretion_radius_cells_.
    ASSERT2("EnzoMethodAccretion",
	    "The diagonal cell width divided by the minimum "
	    "cell width is %g, and accretion_radius_cells "
	    "is %g. The former must be less than twice the "
	    "latter. It is advised to adjust this either "
	    "by changing the dimensions of the root level "
	    "mesh, via the Mesh:root_size parameter, or by "
	    "changing the accretion:accretion_radius_cells "
	    "parameter",
	    diagonal_over_minimum, accretion_radius_cells_,
	    diagonal_over_minimum < 2.0 * accretion_radius_cells_);

    // Check sink particle attributes
    cello::particle_descr()->check_particle_attribute("sink","mass");
    cello::particle_descr()->check_particle_attribute("sink","x");
    cello::particle_descr()->check_particle_attribute("sink","y");
    cello::particle_descr()->check_particle_attribute("sink","z");
    cello::particle_descr()->check_particle_attribute("sink","vx");
    cello::particle_descr()->check_particle_attribute("sink","vy");
    cello::particle_descr()->check_particle_attribute("sink","vz");
    cello::particle_descr()->check_particle_attribute("sink","accretion_rate");

    return;
}
