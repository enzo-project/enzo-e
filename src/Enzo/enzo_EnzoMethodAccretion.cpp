/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretion.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   24 February 2022
/// @brief  Implementation of EnzoMethodAccretion, a base class
///         for "accretion compute" methods. These methods compute
///         the accretion rate onto sink particles, and change the properties
///         of the particles accordingly. Gas density is reduced by setting
///         negative values for the "density_accreted" field. The
///         "accretion_remove_gas" method then subtracts off density_accreted
///         from the gas density field


#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodAccretion::EnzoMethodAccretion
(double accretion_radius_cells,
 double density_threshold,
 double max_mass_fraction,
 bool   conserve_angular_momentum,
 double ang_mom_threshold_radius_cells)
  : Method(),
    accretion_radius_cells_(accretion_radius_cells),
    density_threshold_(density_threshold),
    max_mass_fraction_(max_mass_fraction),
    conserve_angular_momentum_(conserve_angular_momentum),
    ang_mom_threshold_radius_cells_(ang_mom_threshold_radius_cells),
    ir_accretion_(-1)
{
  // Check if density threshold is at least as large as the density floor
  // set by the VL+CT method
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion",
	 "Density threshold must be at least as large as the density "
	 "floor set by the VL+CT method",
	 density_threshold_ >= enzo::config()->method_vlct_density_floor);

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
	  max_mass_fraction_>= 0.0 && max_mass_fraction <= 1.0);

  // Check that ang_mom_threshold_radius_cells_ is between 0.0 and 0.5
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion()",
	 "Method:accretion:ang_mom_threshold_radius_cells must be less than 0.5.",
	 ang_mom_threshold_radius_cells_ > 0.0 && ang_mom_threshold_radius_cells_ < 0.5);

  const int * ghost_depth = enzo::config()->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  // The number of cell widths within the accretion radius cannot be larger than
  // any of the ghost depths
  ASSERT("EnzoMethodAccretion::EnzoMethodAccretion() ",
	 "The accretion radius must be no greater than the ghost depth"
	 "(4 cells by default)",
	 accretion_radius_cells_ <= min_ghost_depth);

  // define required fields
  cello::define_field("density_source");
  cello::define_field("density_source_accumulate");
  cello::define_field("mom_dens_x_source");
  cello::define_field("mom_dens_x_source_accumulate");
  cello::define_field("mom_dens_y_source");
  cello::define_field("mom_dens_y_source_accumulate");
  cello::define_field("mom_dens_z_source");
  cello::define_field("mom_dens_z_source_accumulate");
  cello::define_field("te_dens_source");
  cello::define_field("te_dens_source_accumulate");

  // if sink particles have a "metal_fraction" attribute, then
  // define a "metal_density" field
  int it = cello::particle_descr()->type_index("sink");
  if (cello::particle_descr()->has_attribute(it,"metal_fraction")){
      cello::define_field("metal_density");
      cello::define_field_in_group("metal_density","color");
  }

  // Initial refresh: refresh all fields and sink particles
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
  if (conserve_angular_momentum_){
    refresh_accretion->add_field_src_dst
    ("mom_dens_x_source","mom_dens_x_source_accumulate");
    refresh_accretion->add_field_src_dst
    ("mom_dens_y_source","mom_dens_y_source_accumulate");
    refresh_accretion->add_field_src_dst
    ("mom_dens_z_source","mom_dens_z_source_accumulate");
    refresh_accretion->add_field_src_dst
    ("te_dens_source","te_dens_source_accumulate");
  }

  refresh_accretion->set_callback(CkIndex_EnzoBlock::p_method_accretion_end());

}

void EnzoMethodAccretion::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | accretion_radius_cells_;
  p | density_threshold_;
  p | max_mass_fraction_;
  p | conserve_angular_momentum_;
  p | ang_mom_threshold_radius_cells_;
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
  method->add_source_fields(this);
  compute_done();
  return;
}

//--------------------------------------------------------------------------------------------
void EnzoMethodAccretion::add_source_fields(EnzoBlock * enzo_block) throw()
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
  enzo_float * te_dens_source = (enzo_float*) field.values("te_dens_source");
  enzo_float * te_dens_source_accumulate =
    (enzo_float*) field.values("te_dens_source_accumulate");

  // Loop over all cells, adding on the "source" and "source_accumulate" fields
  // and rescaling where appropriate. "source" and "source_accumulate" fields
  // are then set to zero, ready for the next cycle.
  for (int i = 0; i < mx * my * mz; i++){

    // Update density
    const enzo_float old_dens = density[i];
    density[i] += density_source[i] + density_source_accumulate[i];
    density_source[i] = 0.0;
    density_source_accumulate[i] = 0.0;

    // Update color fields
    EnzoMethodStarMaker::rescale_densities(enzo_block,i,density[i] / old_dens);

    // Update velocity and total energy fields
    const enzo_float new_mom_dens_x =
      old_dens * vx[i] + mom_dens_x_source[i] + mom_dens_x_source_accumulate[i];
    vx[i] = new_mom_dens_x / density[i];
    mom_dens_x_source[i] = 0.0;
    mom_dens_x_source_accumulate[i] = 0.0;
    const enzo_float new_mom_dens_y =
      old_dens * vx[i] + mom_dens_y_source[i] + mom_dens_y_source_accumulate[i];
    vx[i] = new_mom_dens_y / density[i];
    mom_dens_y_source[i] = 0.0;
    mom_dens_y_source_accumulate[i] = 0.0;
    const enzo_float new_mom_dens_z =
	old_dens * vx[i] + mom_dens_z_source[i] + mom_dens_z_source_accumulate[i];
    vx[i] = new_mom_dens_z / density[i];
    mom_dens_z_source[i] = 0.0;
    mom_dens_z_source_accumulate[i] = 0.0;
    const enzo_float new_te_dens =
      old_dens * specific_te[i] + te_dens_source[i] + te_dens_source_accumulate[i];
    specific_te[i] = new_te_dens / density[i];
    te_dens_source[i] = 0.0;
    te_dens_source_accumulate[i] = 0.0;

  } // Loop over cells

  return;
}

// Required
double EnzoMethodAccretion::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodAccretion::do_checks_() throw()
{
    // Check if merge_sinks method precedes accretion_compute method
    ASSERT("EnzoMethodAccretion",
	   "merge_sinks must precede accretion_compute",
	   enzo::problem()->method_precedes("merge_sinks",
					    "accretion_compute"));

    // Check if merging radius is at least twice that of the accretion
    // radius
    ASSERT("EnzoMethodAccretion::EnzoMethodAccretion() ",
	   "Merging radius (Method:merge_sinks:merging_radius_cells "
	   "must be at least twice the accretion radius "
	   "(Method:accretion_compute:accretion_radius).",
	   enzo::config()->method_merge_sinks_merging_radius_cells >=
	   2.0 * accretion_radius_cells_);

    // Check if VL+CT method is being used.
    ASSERT("EnzoMethodAccretion::EnzoMethodAccretion() ",
	   "accretion_compute requires vlct method. (note: at the "
	   "moment this means that cosmology can't be used in "
	   "combination with accretion.",
	   enzo::problem()->method_exists("vlct"));

    return;
}
