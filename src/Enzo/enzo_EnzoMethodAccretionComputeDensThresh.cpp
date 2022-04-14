// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionCompute.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       10 March 2022
/// @brief      Implementation of EnzoMethodAccretionComputeDensThresh, a class
///             from EnzoMethodAccretionCompute.
///             This method reduces the gas density in the accretion zone around
///             a star particle to a value set by density_threshold_,
///             and adds mass and momentum lost by the gas to the star particle.

#include "cello.hpp"
#include "enzo.hpp"

//------------------------------------------------------------------

EnzoMethodAccretionComputeDensThresh::EnzoMethodAccretionComputeDensThresh
(double accretion_radius_cells,
 double density_threshold)
  : EnzoMethodAccretionCompute(accretion_radius_cells,
			       density_threshold)
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
  // which forces blocks containing accreting star particles to be
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
  const int id   = field.field_id("density");
  const int ida  = field.field_id("density_accreted");
  const int idvx  = field.field_id("velocity_x");
  const int idvy  = field.field_id("velocity_y");
  const int idvz  = field.field_id("velocity_z");

  enzo_float * density = (enzo_float*) field.values(id);
  enzo_float * density_accreted = (enzo_float*) field.values(ida);
  enzo_float * vx = (enzo_float*) field.values(idvx);
  enzo_float * vy = (enzo_float*) field.values(idvy);
  enzo_float * vz = (enzo_float*) field.values(idvz);
  
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

  // Get the coordinates of the center of the cell with index 0 (outermost ghost zone cell)
  double xm, ym, zm;
  block->data()->lower(&xm,&ym,&zm);
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  double first_cell_center_x = (0.5 - gx) * hx + xm;
  double first_cell_center_y = (0.5 - gy) * hy + ym;
  double first_cell_center_z = (0.5 - gz) * hz + zm;

  Particle particle = block->data()->particle();
  int it = particle.type_index("star");
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
	  ceil((px[ip*dp] - first_cell_center_x - accretion_radius) / hx);
	const int min_ind_y =
	  ceil((py[ip*dp] - first_cell_center_y - accretion_radius) / hy);
	const int min_ind_z =
	  ceil((pz[ip*dp] - first_cell_center_z - accretion_radius) / hz);
	const int max_ind_x =
	  floor((px[ip*dp] - first_cell_center_x + accretion_radius) / hx);
	const int max_ind_y =
	  floor((py[ip*dp] - first_cell_center_y + accretion_radius) / hy);
	const int max_ind_z =
	  floor((pz[ip*dp] - first_cell_center_z + accretion_radius) / hz);

	// Loop over all cells in this region
	for (int kk = min_ind_z; kk <= max_ind_z; kk++){
	  for (int jj = min_ind_y; jj <= max_ind_y; jj++){
	    for (int ii = min_ind_x; ii <= max_ind_x; ii++){

	      // Check if center of cell is within accretion radius of particle
	      const double cell_center_x = first_cell_center_x + ii * hx;
	      const double cell_center_y = first_cell_center_y + jj * hy;
	      const double cell_center_z = first_cell_center_z + kk * hz;

	      const double r2 =
		  (px[ip*dp] - cell_center_x) * (px[ip*dp] - cell_center_x)
		+ (py[ip*dp] - cell_center_y) * (py[ip*dp] - cell_center_y)
		+ (pz[ip*dp] - cell_center_z) * (pz[ip*dp] - cell_center_z);

	      if (r2 < accretion_radius * accretion_radius){
		const int index = INDEX(ii,jj,kk,mx,my);
		if (density[index] > density_threshold_){
		  
		  // Set density_accreted equal to the density removed
		  // Mass and momentum removed from gas is added to the star particle
		  const enzo_float density_removed = density[index] - density_threshold_;
		  density_accreted[index] = density_removed;
		  
		  const double mass_removed = density_removed * cell_volume;
		  const double momentum_removed_x = mass_removed * vx[index];
		  const double momentum_removed_y = mass_removed * vy[index];
		  const double momentum_removed_z = mass_removed * vz[index];

		  const double old_particle_momentum_x = pmass[ip*dm] * pvx[ip*dv];
		  const double old_particle_momentum_y = pmass[ip*dm] * pvy[ip*dv];
		  const double old_particle_momentum_z = pmass[ip*dm] * pvz[ip*dv];

		  const double new_particle_momentum_x =
		    old_particle_momentum_x + momentum_removed_x;
		  const double new_particle_momentum_y =
		    old_particle_momentum_y + momentum_removed_y;
		  const double new_particle_momentum_z =
		    old_particle_momentum_z + momentum_removed_z;

		  pmass[ip*dm] += mass_removed;
		  pvx[ip*dv] = new_particle_momentum_x / pmass[ip*dm];
		  pvy[ip*dv] = new_particle_momentum_y / pmass[ip*dm];
		  pvz[ip*dv] = new_particle_momentum_z / pmass[ip*dm];
	    
		} // if density[index] > density_threshold				
	      } // if (r2 < accretion_radius * accretion_radius)      
	    }
	  }
	} // Loop over cells in region bounding accretion zone
      } // Loop over particles in this batch
    } // Loop over batches
  } // if (num_particles > 0)
  return;
}
