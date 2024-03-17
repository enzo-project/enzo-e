// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     7 June 2022
/// @brief    Initializer for the Shu Collapse problem as described
///           in Federrath et al 2010, ApJ, 713, 269.

#include "Enzo/initial/initial.hpp"
#include "Enzo/enzo.hpp"
#include "Cello/cello.hpp"

EnzoInitialShuCollapse::EnzoInitialShuCollapse
(int cycle, double time,
 const double center[3],
 const double drift_velocity[3],
 double truncation_radius,
 double nominal_sound_speed,
 double instability_parameter,
 double external_density,
 bool central_sink_exists,
 double central_sink_mass) throw()
  :Initial(cycle, time),
   truncation_radius_(truncation_radius),
   nominal_sound_speed_(nominal_sound_speed),
   instability_parameter_(instability_parameter),
   external_density_(external_density),
   central_sink_exists_(central_sink_exists),
   central_sink_mass_(central_sink_mass)
{
  center_[0] = center[0];
  center_[1] = center[1];
  center_[2] = center[2];

  drift_velocity_[0] = drift_velocity[0];
  drift_velocity_[1] = drift_velocity[1];
  drift_velocity_[2] = drift_velocity[2];

}

void EnzoInitialShuCollapse::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,center_,3);
  PUParray(p,drift_velocity_,3);
  p | truncation_radius_;
  p | nominal_sound_speed_;
  p | instability_parameter_;
  p | external_density_;
  p | central_sink_exists_;
  p | central_sink_mass_;

}

void EnzoInitialShuCollapse::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()

{
  // Check if we have periodic boundary conditions
  int px,py,pz;
  hierarchy->get_periodicity(&px,&py,&pz);
  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	 "Shu Collapse must have periodic boundary conditions.",px && py && pz);

  // Check if the truncation radius is less than half the domain size
  double dxm,dxp,dym,dyp,dzm,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);
  double domain_width_x = dxp - dxm;
  double domain_width_y = dyp - dym;
  double domain_width_z = dzp - dzm;

  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	 "Truncation radius must no more than a quarter of the domain width in all "
	 "dimensions.",
	 (truncation_radius_ <= 0.25 * domain_width_x) &&
	 (truncation_radius_ <= 0.25 * domain_width_y) &&
	 (truncation_radius_ <= 0.25 * domain_width_z));

  // Check that gamma is sufficiently close to unity
  const double gamma = enzo::fluid_props()->gamma();
  ASSERT("EnzoInitialShuCollapse",
	 "The adiabatic index (Field:gamma) must be between 1 and 1.000001.",
	 gamma > 1.0 && gamma < 1.000001);

  // Check that truncation_radius_ is positive
  ASSERT("EnzoInitialShuCollapse",
	 "Initial:shu_collapse:truncation_radius must be larger than 0.0.",
	 truncation_radius_ > 0.0);

  // Check that nominal_sound_speed_ is positive
  ASSERT("EnzoInitialShuCollapse",
	 "Initial:shu_collapse:nominal_sound_speed must be larger than 0.0.",
	 nominal_sound_speed_ > 0.0);

  // Check that instability_parameter_ is greater than 2.0. This ensures that the system
  // actually collapses.
  ASSERT("EnzoInitialShuCollapse",
	 "Initial:shu_collapse:instability_parameter must be greater than 2.0.",
	 instability_parameter_ > 2.0);

  // Check that central_sink_mass_ is non-negative
  ASSERT("EnzoInitialShuCollapse",
	 "Initial:shu_collapse:central_sink_mass must be non-negative.",
	 truncation_radius_ > 0.0);

  // Check if sink_maker method is being used
  ASSERT("EnzoInitialShuCollapse",
	   "If shu_collapse initializer is used, the sink_maker "
	   "method is required.",
	   enzo::problem()->method_exists("sink_maker"));

  // Check if mhd_vlct method is being used
  ASSERT("EnzoInitialShuCollapse",
	 "If shu_collapse initializer is used, the mhd_vlct method is "
	 "required.",
	 enzo::problem()->method_exists("mhd_vlct"));

  // Check that mhd_choice parameter is set to "no_bfield"
  ASSERT("EnzoInitialShuCollapse",
	 "Method:mhd_vlct:mhd_choice must be set to no_bfield",
	 enzo::config()->method_vlct_mhd_choice == "no_bfield");

  // Check that riemann_solver parameter is set to "hllc"
  ASSERT("EnzoInitialShuCollapse",
	 "Method:mhd_vlct:mhd_choice must be set to hllc",
	 enzo::config()->method_vlct_riemann_solver == "hllc");

  if (!block->is_leaf()) return;
  ASSERT("EnzoInitialShuCollapse",
	 "Block does not exist",
	 block != NULL);

  double folded_center_position[3];
  hierarchy->get_folded_position(center_, folded_center_position);

  center_[0] = folded_center_position[0];
  center_[1] = folded_center_position[1];
  center_[2] = folded_center_position[2];

  Field field = block->data()->field();

  // Get Field parameters
  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
  const int m = mx * my * mz;
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);

  // Get cell widths
  double hx,hy,hz;
  block->cell_width(&hx, &hy, &hz);

  // Get pointers to fields
  // It's likely that not all these fields are required for the problem, and that not all of
  // the required fields need to be initialized, but let's just do it for safety.
  enzo_float *  d           = (enzo_float *) field.values ("density");
  enzo_float * dt           = (enzo_float *) field.values ("density_total");
  enzo_float * dp           = (enzo_float *) field.values ("density_particle");
  enzo_float * dpa          = (enzo_float *) field.values ("density_particle_accumulate");
  enzo_float * dg           = (enzo_float *) field.values ("density_gas");
  enzo_float * pm           = (enzo_float *) field.values ("particle_mass");
  enzo_float * po           = (enzo_float *) field.values ("potential");
  enzo_float * po_temp      = (enzo_float *) field.values ("potential_temp");
  enzo_float * po_copy      = (enzo_float *) field.values ("potential_copy");
  enzo_float * specific_te  = (enzo_float *) field.values ("total_energy");
  enzo_float * pressure     = (enzo_float *) field.values ("pressure");
  enzo_float * ax           = (enzo_float *) field.values ("acceleration_x");
  enzo_float * ay           = (enzo_float *) field.values ("acceleration_y");
  enzo_float * az           = (enzo_float *) field.values ("acceleration_z");
  enzo_float * vx           = (enzo_float *) field.values ("velocity_x");
  enzo_float * vy           = (enzo_float *) field.values ("velocity_y");
  enzo_float * vz           = (enzo_float *) field.values ("velocity_z");
  enzo_float * x            = (enzo_float *) field.values ("X");
  enzo_float * x_copy       = (enzo_float *) field.values ("X_copy");
  enzo_float * b            = (enzo_float *) field.values ("B");
  enzo_float * b_copy       = (enzo_float *) field.values ("B_copy");
  enzo_float * ds           = (enzo_float *) field.values ("density_source");
  enzo_float * dsa          = (enzo_float *) field.values ("density_source_accumulate");
  enzo_float * mdxs         = (enzo_float *) field.values ("mom_dens_x_source");
  enzo_float * mdxsa        = (enzo_float *) field.values ("mom_dens_x_source_accumulate");
  enzo_float * mdys         = (enzo_float *) field.values ("mom_dens_y_source");
  enzo_float * mdysa        = (enzo_float *) field.values ("mom_dens_y_source_accumulate");
  enzo_float * mdzs         = (enzo_float *) field.values ("mom_dens_z_source");
  enzo_float * mdzsa        = (enzo_float *) field.values ("mom_dens_z_source_accumulate");

  // For most of the fields, we initialise their values to zero everywhere
  std::fill_n(dt,m,0.0);
  std::fill_n(dp,m,0.0);
  std::fill_n(dpa,m,0.0);
  std::fill_n(dg,m,0.0);
  std::fill_n(pm,m,0.0);
  std::fill_n(po,m,0.0);
  std::fill_n(po_temp,m,0.0);
  std::fill_n(po_copy,m,0.0);
  std::fill_n(pressure,m,0.0);
  std::fill_n(ax,m,0.0);
  std::fill_n(ay,m,0.0);
  std::fill_n(az,m,0.0);
  std::fill_n(x,m,0.0);
  std::fill_n(x_copy,m,0.0);
  std::fill_n(b,m,0.0);
  std::fill_n(b_copy,m,0.0);
  std::fill_n(ds,m,0.0);
  std::fill_n(dsa,m,0.0);
  std::fill_n(mdxs,m,0.0);
  std::fill_n(mdxsa,m,0.0);
  std::fill_n(mdys,m,0.0);
  std::fill_n(mdysa,m,0.0);
  std::fill_n(mdzs,m,0.0);
  std::fill_n(mdzsa,m,0.0);

  // Set velocity
  std::fill_n(vx,m,drift_velocity_[0]);
  std::fill_n(vy,m,drift_velocity_[1]);
  std::fill_n(vz,m,drift_velocity_[2]);

  // Set specific total energy by computing specific kinetic and internal energies
  const enzo_float specific_ke = 0.5 * (drift_velocity_[0] * drift_velocity_[0] +
					drift_velocity_[1] * drift_velocity_[1] +
					drift_velocity_[2] * drift_velocity_[2]);

  // Set initial pressure to be density times square of `nominal_sound_speed_`
  // (i.e., pretend gas is isothermal). This ensures that it has the same (1/r^2) profile as
  // the density field. This means that the specific internal energy is `nominal_sound_speed_`
  // squared over (`gamma` - 1), where `gamma` is slightly bigger than 1. Note that the true
  // sound speed will be `gamma` multiplied by `nominal_sound_speed_` everywhere initially,
  // and will be non-uniform as the system evolves.
  const enzo_float specific_ie = nominal_sound_speed_ * nominal_sound_speed_ / (gamma - 1.0);
  std::fill_n(specific_te,m,specific_ke + specific_ie);

  // Now to initialise the density field
  const double const_G  = enzo::grav_constant_cgs() * enzo::units()->density() *
    enzo::units()->time() * enzo::units()->time();

  const double density_profile_factor =
    instability_parameter_ * nominal_sound_speed_ * nominal_sound_speed_ /
    (4.0 * cello::pi * const_G);

  for (int iz = 0; iz < mz ; iz++){
    const double z = bzm + (iz - gz + 0.5)*hz;
    for (int iy = 0; iy < my; iy++){
      const double y = bym + (iy - gy + 0.5)*hy;
      for (int ix = 0; ix < mx; ix++){
	const double x = bxm + (ix - gx + 0.5)*hx;
	double cell_pos[3] = {x,y,z};
	double npi[3];
	hierarchy->get_nearest_periodic_image(cell_pos,center_,npi);
	const double r2 =
	  (npi[0] - center_[0]) * (npi[0] - center_[0]) +
	  (npi[1] - center_[1]) * (npi[1] - center_[1]) +
	  (npi[2] - center_[2]) * (npi[2] - center_[2]);
	const int i = INDEX(ix,iy,iz,mx,my);
	d[i] =
	  (r2 < truncation_radius_ * truncation_radius_) ?
	  density_profile_factor / r2 :
	  external_density_;
      } //ix
    } //iy
  } //iz

  // If central_particle is true and collapse center is in this block, we
  // add a particle at the collapse center
  if (central_sink_exists_ &&
      block->check_position_in_block(center_[0],center_[1],center_[2]))
    {

      Particle particle = block->data()->particle();

      // Attribute indices
      const int it   = particle.type_index("sink");
      const int ia_m = particle.attribute_index (it, "mass");
      const int ia_x = particle.attribute_index (it, "x");
      const int ia_y = particle.attribute_index (it, "y");
      const int ia_z = particle.attribute_index (it, "z");
      const int ia_vx = particle.attribute_index (it, "vx");
      const int ia_vy = particle.attribute_index (it, "vy");
      const int ia_vz = particle.attribute_index (it, "vz");
      const int ia_copy = particle.attribute_index (it, "is_copy");

      // Attribrute stride lengths
      const int dm   = particle.stride(it, ia_m);
      const int dp   = particle.stride(it, ia_x);
      const int dv   = particle.stride(it, ia_vx);
      const int dcopy = particle.stride(it, ia_copy);

      /// Initialise pointers for particle attribute arrays
      enzo_float * pmass = 0;
      enzo_float * px   = 0;
      enzo_float * py   = 0;
      enzo_float * pz   = 0;
      enzo_float * pvx  = 0;
      enzo_float * pvy  = 0;
      enzo_float * pvz  = 0;
      int64_t * is_copy = 0;

      // ip_block is the index of the particle in the block
      // ip_batch is the index if the particle in its batch
      // ibatch is the index of the batch
      int ip_block = particle.insert_particles(it, 1);
      int ip_batch;
      int ibatch;
      particle.index(ip_block, &ibatch, &ip_batch);
      enzo::simulation()->data_insert_particles(1);

      // Get pointers to particle attribute arrays
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ibatch);
      px    = (enzo_float *) particle.attribute_array(it, ia_x, ibatch);
      py    = (enzo_float *) particle.attribute_array(it, ia_y, ibatch);
      pz    = (enzo_float *) particle.attribute_array(it, ia_z, ibatch);
      pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ibatch);
      pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ibatch);
      pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ibatch);
      is_copy   = (int64_t *) particle.attribute_array(it, ia_copy, ibatch);

      // Now assign values to attributes
      pmass[ip_batch*dm] = central_sink_mass_;
      px[ip_batch*dp] = center_[0];
      py[ip_batch*dp] = center_[1];
      pz[ip_batch*dp] = center_[2];
      pvx[ip_batch*dp] = drift_velocity_[0];
      pvy[ip_batch*dp] = drift_velocity_[1];
      pvz[ip_batch*dp] = drift_velocity_[2];
      is_copy[ip_batch*dcopy] = 0;

    } // Is there are central particle to place in this block?
  return;
}
