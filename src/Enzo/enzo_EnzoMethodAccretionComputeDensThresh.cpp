// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionCompute.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       10 March 2022
/// @brief      Implementation of EnzoMethodAccretionComputeDensThresh, a class
///             from EnzoMethodAccretionCompute.
///             This method reduces the gas density in the accretion zone around
///             a sink particle to
///             max(density_threshold_,(1-max_mass_fraction)*density),
///             and adds mass and momentum lost by the gas to the sink particle.

#include "cello.hpp"
#include "enzo.hpp"

//------------------------------------------------------------------

EnzoMethodAccretionComputeDensThresh::EnzoMethodAccretionComputeDensThresh
(double accretion_radius_cells,
 double density_threshold,
 double max_mass_fraction,
 bool   conserve_angular_momentum)
  : EnzoMethodAccretionCompute(accretion_radius_cells,
			       density_threshold,
			       max_mass_fraction_,
			       conserve_angular_momentum)
{

}

//-------------------------------------------------------------------

void EnzoMethodAccretionComputeDensThresh::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodAccretionCompute::pup(p); // call parent class pup

  return;
}

//--------------------------------------------------------------------

void EnzoMethodAccretionComputeDensThresh::compute (Block * block) throw()
{

  if (enzo::simulation()->cycle() == enzo::config()->initial_cycle)
    do_checks_();

  // Only call compute_ if block is at highest refinement level.
  // Currently this method can only be used if refinement is turned off
  // (unigrid mode), but in future, we will have a refinement condition
  // which forces blocks containing accreting sink particles to be
  // at the highest refinement level, and using this method will
  // require this refinement condition to be activated.
  if (block->level() == enzo::config()->mesh_max_level) {
    this->compute_(block);
  }

  block->compute_done();

  return;
}

//-------------------------------------------------------------------------

void EnzoMethodAccretionComputeDensThresh::compute_(Block * block)

{

  // Get pointers to field data
  Field field = block->data()->field();

  enzo_float * density     = (enzo_float*) field.values("density");
  enzo_float * density_ch  = (enzo_float*) field.values("density_ch");
  enzo_float * vx          = (enzo_float*) field.values("velocity_x");
  enzo_float * vy          = (enzo_float*) field.values("velocity_y");
  enzo_float * vz          = (enzo_float*) field.values("velocity_z");


  enzo_float * vx_ch = conserve_angular_momentum_ ?
    (enzo_float*) field.values("velocity_x_change") : nullptr;
  enzo_float * vy_ch = conserve_angular_momentum_ ?
    (enzo_float*) field.values("velocity_y_change") : nullptr;
  enzo_float * vz_ch = conserve_angular_momentum_ ?
    (enzo_float*) field.values("velocity_z_change") : nullptr;
  enzo_float * ie    = conserve_angular_momentum_ ?
    (enzo_float*) field.values("internal_energy") : nullptr;
  enzo_float * te    = conserve_angular_momentum_ ?
    (enzo_float*) field.values("total_energy") : nullptr;
  enzo_float * te_ch = conserve_angular_momentum_?
    (enzo_float*) field.values("total_energy_change") : nullptr;


  // Get field dimensions
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  // Get cell widths and cell volume
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  const double cell_volume = hx * hy * hz;

  // Set accretion radius to be accretion_radius_cells_ multipled by minimum cell width
  const double min_cell_width = std::min(hx,std::min(hy,hz));
  const double accretion_radius = accretion_radius_cells_ * min_cell_width;

  double xm, ym, zm;
  block->data()->lower(&xm,&ym,&zm);
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  Particle particle = block->data()->particle();
  int it = particle.type_index("sink");
  int num_particles = particle.num_particles(it);

  if (num_particles > 0) {

    // Declare pointers to particle attributes
    enzo_float *pmass, *px, *py, *pz, *pvx, *pvy, *pvz;

    // Get attribute indices
    const int ia_m   = particle.attribute_index (it, "mass");
    const int ia_x   = particle.attribute_index (it, "x");
    const int ia_y   = particle.attribute_index (it, "y");
    const int ia_z   = particle.attribute_index (it, "z");
    const int ia_vx  = particle.attribute_index (it, "vx");
    const int ia_vy  = particle.attribute_index (it, "vy");
    const int ia_vz  = particle.attribute_index (it, "vz");

    // Attribrute stride lengths
    const int dm   = particle.stride(it, ia_m);
    const int dp   = particle.stride(it, ia_x);
    const int dv   = particle.stride(it, ia_vx);

    // Loop over batches
    const int nb = particle.num_batches(it);
    for (int ib=0; ib < nb; ib++){

      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
      px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);
      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

      // Loop over particles in this batch
      const int np = particle.num_particles(it,ib);
      for (int ip=0; ip < np; ip++){

	// Find indices which specify a region which bounds the accretion zone
	const int min_ind_x =
	  ceil((px[ip*dp] - xm - accretion_radius) / hx - 0.5) + gx;
	const int min_ind_y =
	  ceil((py[ip*dp] - ym - accretion_radius) / hy - 0.5) + gy;
	const int min_ind_z =
	  ceil((pz[ip*dp] - zm - accretion_radius) / hz - 0.5) + gz;
	const int max_ind_x =
	  floor((px[ip*dp] - xm + accretion_radius) / hx - 0.5) + gx;
	const int max_ind_y =
	  floor((py[ip*dp] - ym + accretion_radius) / hy - 0.5) + gy;
	const int max_ind_z =
	  floor((pz[ip*dp] - zm + accretion_radius) / hz - 0.5) + gz;

	// keep track of total mass change
	double total_mass_change = 0.0;

	// keep track of the centre of mass
	double x_com = pmass[ip*dm] * px[ip*dp];
	double y_com = pmass[ip*dm] * py[ip*dp];
	double z_com = pmass[ip*dm] * pz[ip*dp];

	// keep track of total momentum change of particle
	double total_p_momentum_change_x = 0.0;
	double total_p_momentum_change_y = 0.0;
	double total_p_momentum_change_z = 0.0;

	// Loop over all cells in this region
	for (int k = min_ind_z; k <= max_ind_z; k++){
	  for (int j = min_ind_y; j <= max_ind_y; j++){
	    for (int i = min_ind_x; i <= max_ind_x; i++){

	      // Check if center of cell is within accretion radius of particle
	      const double cell_center_x = xm + (i - gx) * hx;
	      const double cell_center_y = ym + (j - gy) * hy;
	      const double cell_center_z = zm + (k - gz) * hz;

	      const double r2 =
		  (px[ip*dp] - cell_center_x) * (px[ip*dp] - cell_center_x)
		+ (py[ip*dp] - cell_center_y) * (py[ip*dp] - cell_center_y)
		+ (pz[ip*dp] - cell_center_z) * (pz[ip*dp] - cell_center_z);

	      if (r2 < accretion_radius * accretion_radius){
		const int index = INDEX(i,j,k,mx,my);
		if (density[index] > density_threshold_){

		  const enzo_float density_removed =
		    std::min(density[index] - density_threshold_,
			     max_mass_fraction_ * density[index]);

		  // Set density_change equal to (minus) the density removed
		  density_ch[index] = -density_removed;

		  const double mass_removed = density_removed * cell_volume;

		  // variables used to store the momentum change of the particle
		  // due to this cell being accreted
		  double p_momentum_change_x;
		  double p_momentum_change_y;
		  double p_momentum_change_z;
		  if (conserve_angular_momentum_){

		    // Get radial component of velocity
		    // Technically, what I call the "radial velocity"
		    // is really the radial velocity in the frame of reference
		    // where the particle is at rest, plus the particle's
		    // velocity.
		    const double r = sqrt(r2);
		    const double vx_radial =
		      (vx[index] - pvx[ip*dv]) * (cell_center_x - px[ip*dp]) / r
		      + pvx[ip*dv];
		    const double vy_radial =
		      (vy[index] - pvy[ip*dv]) * (cell_center_x - px[ip*dp]) / r
		      + pvy[ip*dv];
		    const double vz_radial =
		      (vz[index] - pvz[ip*dv]) * (cell_center_x - px[ip*dp]) / r
		      + pvz[ip*dv];

		    // Only the radial momentum of gas is changed (due to mass loss).
		    // This momentum is added to particle
		    p_momentum_change_x = mass_removed * vx_radial;
		    p_momentum_change_y = mass_removed * vy_radial;
		    p_momentum_change_z = mass_removed * vz_radial;

		    // gas momentum change is minus the particle momentum change.
		    // Use this to calculate new velocities
		    const double old_gas_mass = density[index] * cell_volume;
		    const double new_gas_mass = old_gas_mass - mass_removed;
		    const double vx_new =
		      (old_gas_mass * vx[index] - p_momentum_change_x) / new_gas_mass;
		    const double vy_new =
		      (old_gas_mass * vy[index] - p_momentum_change_y) / new_gas_mass;
		    const double vz_new =
		      (old_gas_mass * vz[index] - p_momentum_change_z) / new_gas_mass;

		    // Set "velocity change" fields to account for changes in velocities
		    vx_ch[index] = vx_new - vx[index];
		    vy_ch[index] = vy_new - vy[index];
		    vz_ch[index] = vz_new - vz[index];

		    // Also need to account for change in total energy (really specific
		    // energy)
		    // Not sure if this depends on whether "dual energy" is being used.
		    const double new_total_energy =
		      ie[index] + 0.5 * (vx_new * vx_new + vy_new * vy_new + vz_new * vz_new);
		    te_ch[index] = new_total_energy - te[index];


		  } else {

		    // In this case, we don't change the velocity of the gas, so
		    // momentum change of particle just depends on the mass removed
		    // from the gas
		    p_momentum_change_x = mass_removed * vx[index];
		    p_momentum_change_y = mass_removed * vy[index];
		    p_momentum_change_z = mass_removed * vz[index];

		  }

		  // Update total mass change
		  total_mass_change += mass_removed;

		  // Update center of mass
		  x_com += mass_removed * cell_center_x;
		  y_com += mass_removed * cell_center_y;
		  z_com += mass_removed * cell_center_z;

		  // Update change of particle momentum
		  total_p_momentum_change_x += p_momentum_change_x;
		  total_p_momentum_change_x += p_momentum_change_y;
		  total_p_momentum_change_x += p_momentum_change_z;
		}
	      } // if density[index] > density_threshold
	    } // if (r2 < accretion_radius * accretion_radius)
	  }
	} // Loop over cells in region bounding accretion zone 

	// Update particle properties

	const double old_p_momentum_x = pmass[ip*dm] * pvx[ip*dv];
	const double old_p_momentum_y = pmass[ip*dm] * pvy[ip*dv];
	const double old_p_momentum_z = pmass[ip*dm] * pvz[ip*dv];
	const double new_p_momentum_x =
	  old_p_momentum_x + total_p_momentum_change_x;
	const double new_p_momentum_y =
	  old_p_momentum_y + total_p_momentum_change_y;
	const double new_p_momentum_z =
	  old_p_momentum_z + total_p_momentum_change_z;

	pmass[ip*dm] += total_mass_change;
	px[ip*dp] = x_com / pmass[ip*dm];
	py[ip*dp] = y_com / pmass[ip*dm];
	pz[ip*dp] = z_com / pmass[ip*dm];
	pvx[ip*dv] = new_p_momentum_x / pmass[ip*dm];
	pvy[ip*dv] = new_p_momentum_y / pmass[ip*dm];
	pvz[ip*dv] = new_p_momentum_z / pmass[ip*dm];
      } // Loop over particles in this batch
    } // Loop over batches
  } // if (num_particles > 0)
  return;
}
