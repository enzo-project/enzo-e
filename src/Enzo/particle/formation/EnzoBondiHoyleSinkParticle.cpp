/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoBondiHoyleSinkParticle.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @author John Regan (john.regan@mu.ie)
/// @date   17 May 2022
/// @brief  Implementation of EnzoBondiHoyleSinkParticle, a class which derives
///         from EnzoSinkParticle. Contains additional data and methods used for
///         computing accretion rates according to the method described in
///         Krumholz+ 2004, ApJ, 611, 399


#include "cello.hpp"
#include "enzo.hpp"

//-------------------------------------------------------------------------------

EnzoBondiHoyleSinkParticle::EnzoBondiHoyleSinkParticle
(Block * block,
 int ib,
 int ip,
 double accretion_radius,
 double const_G)
  : EnzoSinkParticle(block,ib,ip,accretion_radius),
    const_G_(const_G)
{
  init_();
}

//--------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::init_() throw()
{
  // Get the field dimensions
  Field field = block_->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  set_host_cell_indices_();
  set_v_inf_2_();
  set_c_s_inf_2_();
  set_bondi_hoyle_radius_();

  // Will create a vector of indices of cells in the accretion zone.
  // Need to keep track of which entry in the vector corresponds to
  // the host cell index
  decltype(acc_zone_1d_indices_.size()) current_vector_index = 0;
  decltype(acc_zone_1d_indices_.size()) host_cell_vector_index;

  // Keep track of the sum of (un-normalized) cell weights,
  // will use this to normalize the weights later.
  double sum_of_cell_weights = 0.0;

  // The bound fraction of the host cell will be the minimum
  // of the bound fractions of all 26 neighboring cells.
  double min_neighboring_bound_fraction = 1.0;

  // Loop over cells in region bounding accretion zone
  for (int iz = min_ind_z_; iz < max_ind_z_; iz++){
    for (int iy = min_ind_y_; iy < max_ind_y_; iy++){
      for (int ix = min_ind_x_; ix < max_ind_x_; ix++){
	double r2;
	if (cell_in_accretion_zone(ix,iy,iz,&r2)) {
	  const int index = INDEX(ix,iy,iz,mx,my);

	  // append index to vector of accretion zone cell indices
	  acc_zone_1d_indices_.push_back(index);

	  // Compute cell weight
	  const double cell_weight = get_cell_weight_(r2);
	  cell_weights_.push_back(cell_weight);
	  sum_of_cell_weights += cell_weight;

	  // Now to compute the bound fraction.
	  // If this is the host cell, set equal to zero, and set
	  // host_cell_vector_index equal to current_vector_index
	  // If this is a cell neighboring the host cell, update
	  // min_neighboring_bound_fraction.
	  double bound_fraction = 0.0;
	  if (index == host_cell_1d_index_)
	      host_cell_vector_index = current_vector_index;
	  else {
	    bound_fraction = get_cell_bound_fraction_(ix,iy,iz);

	    // This cell neighbors the host cell if the maximum difference
	    // between the x/y/z index if this cell and the x/y/z index
	    // of the host cell is 1
	    if (std::max(std::max(std::abs(ix - host_cell_ix_),
				  std::abs(iy - host_cell_iy_)),
			   std::abs(iz - host_cell_iz_)))
	      min_neighboring_bound_fraction =
		std::min(min_neighboring_bound_fraction,bound_fraction);
	  }
	  bound_fractions_.push_back(bound_fraction);
	  current_vector_index++;
	} // if cell in accretion zone
      }
    }
  } // Triple loop over cells in region bounding accretion zone

  // Normalize cell weights
  for (auto &w : cell_weights_) w /= sum_of_cell_weights;

  // Set bound fraction of host cell
  bound_fractions_[host_cell_vector_index] = min_neighboring_bound_fraction;

  // Compute accretion rate
  set_accretion_rate_();

  return;
}

//---------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::set_host_cell_indices_() throw()
{
  Field field = block_->data()->field();
  // Get field dimensions
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  // Get cell widths
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);

  // Get the coordinates of the corner of the block
  double xm, ym, zm;
  block_->data()->lower(&xm,&ym,&zm);

  // Get the ghost depths
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  host_cell_ix_ = floor((px_ - xm) / hx) + gx;
  host_cell_iy_ = floor((py_ - ym) / hy) + gy;
  host_cell_iz_ = floor((pz_ - zm) / hz) + gz;

  host_cell_1d_index_ =
    INDEX(host_cell_ix_,host_cell_iy_,host_cell_iz_,mx,my);

  return;
}

//-----------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::set_bondi_hoyle_radius_() throw()
{
  // Compute the Bondi-Hoyle radius (Equation 10 in Krumholz paper)
  r_BH_ = const_G_ * pmass_ / (v_inf_2_ + c_s_inf_2_);

  return;
}

//----------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::set_v_inf_2_() throw()
{
  Field field = block_->data()->field();
  enzo_float * vx = (enzo_float*) field.values("velocity_x");
  enzo_float * vy = (enzo_float*) field.values("velocity_y");
  enzo_float * vz = (enzo_float*) field.values("velocity_z");

  v_inf_2_ =
    (vx[host_cell_1d_index_] - pvx_) * (vx[host_cell_1d_index_] - pvx_) +
    (vy[host_cell_1d_index_] - pvy_) * (vy[host_cell_1d_index_] - pvy_) +
    (vz[host_cell_1d_index_] - pvz_) * (vz[host_cell_1d_index_] - pvz_);

  return;
}

//------------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::set_c_s_inf_2_() throw()
{
  Field field = block_->data()->field();

  // Get pointers to density, velocity and total energy fields
  const enzo_float * density     = (enzo_float*) field.values("density");
  const enzo_float * vx          = (enzo_float*) field.values("velocity_x");
  const enzo_float * vy          = (enzo_float*) field.values("velocity_y");
  const enzo_float * vz          = (enzo_float*) field.values("velocity_z");
  const enzo_float * specific_te = (enzo_float*) field.values("total_energy");

  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();

  // Check if dual-energy formalism is in use. If so, then this needs to use
  // then the simulation evolves a "specific internal energy" field.
  const bool dual_energy = ! fluid_props->dual_energy_config().is_disabled();
  const enzo_float * specific_ie_field
    = dual_energy ? (enzo_float*) field.values("internal_energy") : nullptr;

  // Get pointers to magnetic field values if they exist.
  enzo_float * bx =
    field.is_field("bfield_x") ? (enzo_float*) field.values("bfield_x") : nullptr;
  enzo_float * by =
    field.is_field("bfield_y") ? (enzo_float*) field.values("bfield_y") : nullptr;
  enzo_float * bz =
    field.is_field("bfield_z") ? (enzo_float*) field.values("bfield_z") : nullptr;

  // Copy host_cell_1d_index_ into a new const int called `i` for make code less verbose
  const int i = host_cell_1d_index_;

  // Get the specific internal energy
  enzo_float specific_ie_cell;
  if (dual_energy) specific_ie_cell = specific_ie_field[i];
  else if (bx)
    specific_ie_cell = specific_te[i] -
      0.5 * ( vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] ) -
      0.5 * ( bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i] ) / density[i];
  else specific_ie_cell = specific_te[i] -
	 0.5 * ( vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] );

  // Now compute and return the square of the sound speed
  const double gamma = fluid_props->gamma();
  c_s_inf_2_ =  gamma * (gamma - 1.0) * specific_ie_cell;

  return;
}

//-----------------------------------------------------------------------------------

double EnzoBondiHoyleSinkParticle::get_cell_weight_(double r2) throw()
{
  // Get cell widths
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);

  // Get mininumum cell width
  double min_cell_width = std::min(std::min(hx,hy),hz);

  // Set the kernel radius (Equation 13 of Krumholz paper)
  const double kernel_radius =
    (r_BH_ < 0.25 * min_cell_width) ? 0.25 * min_cell_width :
    ((r_BH_ <= 0.5 * accretion_radius_) ? r_BH_ : 0.5 * accretion_radius_);

  const double kernel_radius_2 = kernel_radius * kernel_radius;

  // Compute weight from kernel radius (Equation 14).
  return exp(-r2 / kernel_radius_2);
}

//----------------------------------------------------------------------------------

double EnzoBondiHoyleSinkParticle::get_cell_bound_fraction_
(int ix, int iy, int iz) throw()
{
  // Get cell widths
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);

  // Get mininumum cell width
  double min_cell_width = std::min(std::min(hx,hy),hz);

  // Get the coordinates of the corner of the block
  double xm, ym, zm;
  block_->data()->lower(&xm,&ym,&zm);

  // Divide cell into 8^3 (512) "mini-cells".
  // Need the coordinates of the center of the first mini-cell (`x0`,`y0`,`z0`),
  // and the spacing between the mini-cells (`dx`, `dy`, `dz`).
  double x0 = xm + (ix + 0.0625) * hx; // 0.0625 = 1/16
  double y0 = ym + (iy + 0.0625) * hy;
  double z0 = zm + (iz + 0.0625) * hz;

  const double dx = 0.125 * hx; // 0.125 = 1/8
  const double dy = 0.125 * hy;
  const double dz = 0.125 * hz;

  // Initialise bound_fraction to 0 and bound_fraction_increment to (1/8)^3
  double bound_fraction = 0.0;
  double bound_fraction_increment = 0.125 * 0.125 * 0.125;

  // Loop over mini-cells.
  for (int i_mc = 0; i_mc < 8; i_mc++){
    for (int j_mc = 0; j_mc < 8; j_mc++){
      for (int k_mc = 0; k_mc < 8; k_mc++){

	// Compute coordinates of center of mini-cell
	const double x = x0 + i_mc * dx;
	const double y = y0 + j_mc * dy;
	const double z = z0 + k_mc * dz;

	// Boolean which stores whether the mini-cell is bound to the sink particle
	bool mini_cell_is_bound;

	// Compute square of distance of center of mini-cell from sink particle
	const double r2 =
	  (x - px_) * (x - px_) + (y - py_) * (y - py_) + (z - pz_) * (z - pz_);

	// If r is less than 10^(-6) times the minimum cell width, assume mini-cell is
	// bounded. This avoids having to compute gravitational potential energy for small
	// r. In practice, this should never happen, since this function is not called for
	// the sink particle's host cell, but let's just put this in for safety.
	if (r2 < 1.0e-12 * min_cell_width * min_cell_width) mini_cell_is_bound = true;

	else {
	  // Treat the mini-cell as a test particle on a Keplerian orbit around the sink
	  // particle.
	  // If its specific orbital energy is non-negative, it is not bound to the particle.
	  // If negative, we compute the periapsis of the orbit, and if it is less than
	  // a quarter of the mininum cell width, then the mini-cell is bound (Equation 15).
	  // Note that if the current position is less than a quarter of the minimum cell
	  // width from the sink particle, and the specific orbital energy is negative,
	  // then we do not need to compute periapsis since we already know it must be less than
	  // or equal the current distance to the sink particle.
	  const double r = sqrt(r2);
	  double vx_rel, vy_rel, vz_rel;
	  compute_relative_velocity_(x,y,z,&vx_rel,&vy_rel,&vz_rel);

	  const double specific_orbital_energy =
	    0.5 * (vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel)
	    - const_G_ * pmass_ / r ;

	  if (specific_orbital_energy >= 0.0) mini_cell_is_bound = false;
	  else if (r <= 0.25 * min_cell_width) mini_cell_is_bound = true;
	  else {

	    // Compute angular momentum (cross product of `r` and `v`)
	    const double j_x = y * vz_rel - z * vy_rel;
	    const double j_y = z * vx_rel - x * vz_rel;
	    const double j_z = x * vy_rel - y * vx_rel;
	    const double j_mag_2 =
	      j_x * j_x + j_y * j_y + j_z * j_z;

	    // Compute periapsis from Equation 15
	    // Note: there is an error in Equation 15 as written in the
	    // paper - the "specific angular momentum" should be the
	    // "square of the specific angular momentum"
	    const double r_min =
	      -const_G_ * pmass_ / specific_orbital_energy *
	      (1.0 - sqrt(1.0 + 2.0 * j_mag_2 * specific_orbital_energy /
			  (const_G_ * const_G_ * pmass_ * pmass_)));

	    mini_cell_is_bound = (r_min <= 0.25 * min_cell_width);

	  }
	}

	if (mini_cell_is_bound) bound_fraction += bound_fraction_increment;
      }
    }
  } // Loop over mini-cells

  return bound_fraction;
}

//-----------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::set_accretion_rate_() throw()
{
  // Get pointer to density field data
  Field field = block_->data()->field();
  enzo_float * density = (enzo_float*) field.values("density");

  // Compute the weighted mean density
  double weighted_mean_density = 0.0;
  for (decltype(cell_weights_.size()) i = 0; i < cell_weights_.size(); ++i)
    weighted_mean_density += cell_weights_[i] * density[acc_zone_1d_indices_[i]];

  // Get min cell width
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  const double min_cell_width = std::min(hx,std::min(hy,hz));

  // Compute rho_inf using weighted mean density and "alpha function"
  // (Equation 12 of Krumholz paper)
  const double rho_inf = weighted_mean_density /
    alpha_(1.2 * min_cell_width / r_BH_);

  // Now compute accretion rate from Equation 11 of Krumholz paper
  // Value of lambda_ assumes gas is isothermal.
  const double lambda_c = 0.25 * exp(1.5);
  accretion_rate_ = 4.0 * cello::pi * rho_inf * r_BH_ * r_BH_ *
    sqrt(lambda_c * lambda_c * c_s_inf_2_ + v_inf_2_);

  return;

}

//-------------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::compute_relative_velocity_
(const double x,const double y,const double z,
 double* vx_rel, double* vy_rel,double* vz_rel) throw()
{
  // Get block extents and cell widths
  double xm, ym, zm;
  double xp, yp, zp;
  double hx, hy, hz;
  block_->lower(&xm, &ym, &zm);
  block_->upper(&xp, &yp, &zp);
  block_->cell_width(&hx, &hy, &hz);

  // Get field dimensions
  Field field = block_->data()->field();
  int mx, my, mz;
  field.dimensions(0, &mx, &my, &mz);
  int nx, ny, nz;
  field.size(&nx, &ny, &nz);
  int gx, gy, gz;
  field.ghost_depth(0, &gx, &gy, &gz);

  // Get pointers to the gas velocity fields
  enzo_float * vx_gas_field  = (enzo_float*) field.values("velocity_x");
  enzo_float * vy_gas_field  = (enzo_float*) field.values("velocity_y");
  enzo_float * vz_gas_field  = (enzo_float*) field.values("velocity_z");

  // Use linear interpolation to compute gas velocity at position (`x`,`y`,`z`).
  // Partially copied from EnzoMethodPmDeposit
  const double tx = nx * (x - xm) / (xp - xm) - 0.5;
  const double ty = ny * (y - ym) / (yp - ym) - 0.5;
  const double tz = nz * (z - zm) / (zp - zm) - 0.5;

  const int ix0 = gx + floor(tx);
  const int iy0 = gy + floor(ty);
  const int iz0 = gz + floor(tz);

  const int ix1 = ix0 + 1;
  const int iy1 = iy0 + 1;
  const int iz1 = iz0 + 1;

  const double x0 = 1.0 - (tx - floor(tx));
  const double y0 = 1.0 - (ty - floor(ty));
  const double z0 = 1.0 - (tz - floor(tz));

  const double x1 = 1.0 - x0;
  const double y1 = 1.0 - y0;
  const double z1 = 1.0 - z0;

  const int i000 = INDEX(ix0,iy0,iz0,mx,my);
  const int i001 = INDEX(ix0,iy0,iz1,mx,my);
  const int i010 = INDEX(ix0,iy1,iz0,mx,my);
  const int i011 = INDEX(ix0,iy1,iz1,mx,my);
  const int i100 = INDEX(ix1,iy0,iz0,mx,my);
  const int i101 = INDEX(ix1,iy0,iz1,mx,my);
  const int i110 = INDEX(ix1,iy1,iz0,mx,my);
  const int i111 = INDEX(ix1,iy1,iz1,mx,my);

  const double vx_gas =
      x0 * y0 * z0 * vx_gas_field[i000] +  x0 * y0 * z1 * vx_gas_field[i001]
    + x0 * y1 * z0 * vx_gas_field[i010] +  x0 * y1 * z1 * vx_gas_field[i011]
    + x1 * y0 * z0 * vx_gas_field[i100] +  x1 * y0 * z1 * vx_gas_field[i101]
    + x1 * y1 * z0 * vx_gas_field[i110] +  x1 * y1 * z1 * vx_gas_field[i111];

  const double vy_gas =
      x0 * y0 * z0 * vy_gas_field[i000] +  x0 * y0 * z1 * vy_gas_field[i001]
    + x0 * y1 * z0 * vy_gas_field[i010] +  x0 * y1 * z1 * vy_gas_field[i011]
    + x1 * y0 * z0 * vy_gas_field[i100] +  x1 * y0 * z1 * vy_gas_field[i101]
    + x1 * y1 * z0 * vy_gas_field[i110] +  x1 * y1 * z1 * vy_gas_field[i111];

  const double vz_gas =
      x0 * y0 * z0 * vz_gas_field[i000] +  x0 * y0 * z1 * vz_gas_field[i001]
    + x0 * y1 * z0 * vz_gas_field[i010] +  x0 * y1 * z1 * vz_gas_field[i011]
    + x1 * y0 * z0 * vz_gas_field[i100] +  x1 * y0 * z1 * vz_gas_field[i101]
    + x1 * y1 * z0 * vz_gas_field[i110] +  x1 * y1 * z1 * vz_gas_field[i111];

  // Set relative velocity to be gas velocity minus sink particle velocity
  *vx_rel = vx_gas - pvx_;
  *vy_rel = vy_gas - pvy_;
  *vz_rel = vz_gas - pvz_;

  return;
}

//-------------------------------------------------------------------------------------

/* Ported (copy and paste) direct from enzo-dev */
double EnzoBondiHoyleSinkParticle::alpha_(double x) throw()
{

#define XMIN 0.01
#define XMAX 2.0
#define NTABLE 51

  double lambda_c, xtable, xtablep1, alpha_exp;
  int idx;

  /* This is a precomputed table of alpha values.  These correspond to x values
     that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
     this choice of range is that the asymptotic expressions are accurate to
     better than 2% outside this range */

  double alphatable[NTABLE] = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,
			      277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613,
			      82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989,
			      29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714,
			      11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812,
			      4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612,
			      2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,
			      1.62439};

  // A constant that appears in the following formulae.  This hardcoded value is
  // valid for an isothermal gas.
  lambda_c = 0.25*exp(1.5);

  // deal with the off-the-table cases
  if (x < XMIN){
    return lambda_c / sqrt(2.*x*x);
  }

  else if (x >= XMAX){
    return exp(1./x);
  }
  else {
    // we are on the table
    idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
    xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
    xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
    alpha_exp = log(x/xtable) / log(xtablep1/xtable);

    return alphatable[idx] * pow(alphatable[idx+1]/alphatable[idx],alpha_exp);
  }

#undef NTABLE
#undef XMIN
#undef XMAX

}

//------------------------------------------------------------------------------------------

void EnzoBondiHoyleSinkParticle::compute(double density_threshold,
					 double max_mass_fraction) throw()
{
  // Check whether vectors have non-zero size and are all the same size
  ASSERT("EnzoBondiHoyleSinkParticle::compute_()",
	 "Function has been called when acc_zone_1d_indices_ has size of zero.",
	 acc_zone_1d_indices_.size() != 0);

  ASSERT("EnzoBondiHoyleSinkParticle::compute_()",
	 "Function has been called when cell_weights_ has size of zero.",
	 cell_weights_.size() != 0);

  ASSERT("EnzoBondiHoyleSinkParticle::compute_()",
	 "Function has been called when bound_fractions_ has size of zero.",
	 bound_fractions_.size() != 0);

  ASSERT3("EnzoBondiHoyleSinkParticle::compute_()",
	  "acc_zone_1d_indices_ has size %d. \n"
	  "cell_weights_ has size %d. \n"
	  "bound_fractions_ has size %d. \n",
	  acc_zone_1d_indices_.size(),
	  cell_weights_.size(),
	  bound_fractions_.size(),
	  (acc_zone_1d_indices_.size() == cell_weights_.size())
	  && (acc_zone_1d_indices_.size() == bound_fractions_.size()));

  const int n_cells = acc_zone_1d_indices_.size();

  // Get cell volume
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  const double inv_cell_volume = 1.0 / (hx * hy * hz);

  // Get pointer to density field data
  Field field = block_->data()->field();
  enzo_float * density = (enzo_float*) field.values("density");

  // Loop over cells in accretion zone
  for (int i = 0; i < n_cells; i++){
    const int index = acc_zone_1d_indices_[i];
    double density_change =
      bound_fractions_[i] * cell_weights_[i] * accretion_rate_ * block_->state()->dt() * inv_cell_volume;

    // Restrict density change via density_threshold and max_mass_fraction
    density_change = std::min(std::min(density_change,density[index] - density_threshold),
			      max_mass_fraction * density[index]);

    update(density_change,index);

  } // Loop over cells in accretion zone

  return;
}
