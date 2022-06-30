// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialBBTest.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     10 June 2022
/// @brief    Initializer for the "BB Test" problem as described
///           in Federrath et al 2010, ApJ, 713, 269.

#include "cello.hpp"
#include "enzo.hpp"

EnzoInitialBBTest::EnzoInitialBBTest
(int cycle, double time,
 const double center[3],
 const double drift_velocity[3],
 double mean_density,
 double fluctuation_amplitude,
 double truncation_radius,
 double nominal_sound_speed,
 double angular_rotation_velocity,
 double external_density) throw()
  :Initial(cycle, time),
   mean_density_(mean_density),
   fluctuation_amplitude_(fluctuation_amplitude),
   truncation_radius_(truncation_radius),
   nominal_sound_speed_(nominal_sound_speed),
   angular_rotation_velocity_(angular_rotation_velocity),
   external_density_(external_density)
{
  center_[0] = center[0];
  center_[1] = center[1];
  center_[2] = center[2];

  drift_velocity_[0] = drift_velocity[0];
  drift_velocity_[1] = drift_velocity[1];
  drift_velocity_[2] = drift_velocity[2];
}

void EnzoInitialBBTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,center_,3);
  PUParray(p,drift_velocity_,3);
  p | mean_density_;
  p | fluctuation_amplitude_;
  p | truncation_radius_;
  p | nominal_sound_speed_;
  p | angular_rotation_velocity_;
  p | external_density_;
}

void EnzoInitialBBTest::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()

{
  // Check if the truncation radius is less than half the domain size
  double dxm,dxp,dym,dyp,dzm,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);
  double domain_width_x = dxp - dxm;
  double domain_width_y = dyp - dym;
  double domain_width_z = dzp - dzm;

  ASSERT("EnzoInitialBBTest::EnzoInitialBBTest()",
	 "Truncation radius must be no more than one quarter of the domain width in "
	 "all dimensions.",
	 (truncation_radius_ <= 0.25 * domain_width_x) &&
	 (truncation_radius_ <= 0.25 * domain_width_y) &&
	 (truncation_radius_ <= 0.25 * domain_width_z));

  // Check that gamma is sufficiently close to unity
  const double gamma = enzo::fluid_props()->gamma();
  ASSERT("EnzoInitialBBTest",
	 "The adiabatic index (Field:gamma) must be between 1 and 1.000001.",
	 gamma > 1.0 && gamma < 1.000001);

  // Check that truncation_radius_ is positive
  ASSERT("EnzoInitialBBTest",
	 "Initial:bb_test:truncation_radius must be larger than 0.0.",
	 truncation_radius_ > 0.0);

  // Check that nominal_sound_speed_ is positive
  ASSERT("EnzoInitialBBTest",
	 "Initial:bb_test:nominal_sound_speed must be larger than 0.0.",
	 nominal_sound_speed_ > 0.0);

  // Check if mhd_vlct method is being used
  ASSERT("EnzoInitialBBTest",
	 "If bb_test initializer is used, the mhd_vlct method is "
	 "required.",
	 enzo::problem()->method_exists("mhd_vlct"));

  // Check that mhd_choice parameter is set to "no_bfield"
  ASSERT("EnzoInitialBBTest",
	 "Method:mhd_vlct:mhd_choice must be set to no_bfield",
	 enzo::config()->method_vlct_mhd_choice == "no_bfield");

  // Check that riemann_solver parameter is set to "hllc"
  ASSERT("EnzoInitialBBTest",
	 "Method:mhd_vlct:mhd_choice must be set to hllc",
	 enzo::config()->method_vlct_riemann_solver == "hllc");

  // Check that (1 - fluctuation_amplitude_) * mean_density is larger
  // than density floor
  //
  // TODO: remove the use of density_dbl_prec. This is a temporary workaround
  ASSERT("EnzoInitialAccretionTest",
	 "Initial gas density must be at least as large as the density "
	 "floor set by the ppm method",
	 (1.0 - fluctuation_amplitude_) * mean_density_ >=
         enzo::fluid_props()->fluid_floor_config().density_dbl_prec());

  // Check if sink_maker method is being used
  ASSERT("EnzoInitialBBTest",
	 "If bb_test initializer is used, the sink_maker "
	 "method is required.",
	 enzo::problem()->method_exists("sink_maker"));

  if (!block->is_leaf()) return;
  ASSERT("EnzoInitialBBTest",
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

  for (int iz = 0; iz < mz ; iz++){
    const double z = bzm + (iz - gz + 0.5)*hz;
    for (int iy = 0; iy < my; iy++){
      const double y = bym + (iy - gy + 0.5)*hy;
      for (int ix = 0; ix < mx; ix++){
	const int i = INDEX(ix,iy,iz,mx,my);
	const double x = bxm + (ix - gx + 0.5)*hx;
	double cell_pos[3] = {x,y,z};
	double npi[3];
	hierarchy->get_nearest_periodic_image(cell_pos,center_,npi);
	const double x_rel = npi[0] - center_[0];
	const double y_rel = npi[1] - center_[1];
	const double z_rel = npi[2] - center_[2];
	const double r2 = x_rel * x_rel + y_rel * y_rel + z_rel * z_rel;

	if (r2 > truncation_radius_ * truncation_radius_) {

	  // If outside truncation radius, set density to external_density_,
	  // and velocity equal to drift velocity.
	  d[i]  = external_density_;
	  vx[i] = drift_velocity_[0];
	  vy[i] = drift_velocity_[1];
	  vz[i] = drift_velocity_[2];

	} else if (r2 < 1.0e-12 * hx * hx) {
	  d[i]  = mean_density_;
	  vx[i] = drift_velocity_[0];
	  vy[i] = drift_velocity_[1];
	  vz[i] = drift_velocity_[2];

	} else {

	  double tiny_number = 1.0e-100 * hx;

	  // Calculate spherical coordinates `r`, `theta` and `phi`
	  if (std::abs(x_rel < tiny_number) && (y_rel == 0.0)){
	    // In this case, where x_rel is tiny and y_rel is zero, phi is undefined,
	    // so just set density to be mean density, set velocity to be drift velocity
	    d[i]  = mean_density_;
	    vx[i] = drift_velocity_[0];
	    vy[i] = drift_velocity_[1];
	    vz[i] = drift_velocity_[2];
	  } else {

	    const double r = sqrt(r2);
	    const double theta = acos(z_rel / r);
	    double phi;

	    // computation of phi depends on which quadrant we are in
	    if (x_rel > tiny_number) phi = atan(y/x);
	    else if (x_rel < -tiny_number) {
	      phi = (y_rel >= 0.0) ? (atan(y/x) + cello::pi) : (atan(y/x) - cello::pi);
	    } else {
	      // Magnitude of x_rel is less than tiny number
	      // Case where y_rel is zero has already been handled;
	      phi = (y_rel > 0.0) ? (0.5 * cello::pi) : (-0.5 * cello::pi);
	    }

	    // density is determined by phi
	    d[i] = mean_density_ * (1.0 + fluctuation_amplitude_ * cos(2.0 * phi));

	    // solid body rotation around the z-axis
	    // z-component of velocity is always just the drift velocity
	    vz[i] = drift_velocity_[2];

	    // x-component
	    vx[i] = -angular_rotation_velocity_ * r * sin(theta) * sin(phi) + drift_velocity_[0];

	    // y-component
	    vy[i] =  angular_rotation_velocity_ * r * sin(theta) * cos(phi) + drift_velocity_[1];

	  } // Close else on line 253
	} // Close else on line 241

	// Set specific total energy by computing specific kinetic and internal energies
	const enzo_float specific_ke =
	  0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

	// Set initial pressure to be density times square of `nominal_sound_speed_`
	// (i.e., pretend gas is isothermal). This ensures that it has the same (1/r^2) profile
	// as the density field. This means that the specific internal energy is
	// `nominal_sound_speed_` squared over (`gamma` - 1), where `gamma` is slightly bigger
	// than 1. Note that the true sound speed will be `gamma` multiplied by
	// `nominal_sound_speed_` everywhere initially, and will be non-uniform as the
	// system evolves.
	const enzo_float specific_ie =
	  nominal_sound_speed_ * nominal_sound_speed_ / (gamma - 1.0);

	specific_te[i] = specific_ke + specific_ie;
      } //ix
    } //iy
  } //iz
  return;
}
