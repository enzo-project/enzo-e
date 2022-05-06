/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretion.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   5 May 2022
/// @brief  Implementation of EnzoSinkParticle, a class which encapsulates
///         the data associated with a sink particle and its accretion zone,
///         as well as associated methods for computing accretion rates and
///         reading and writing to / from the particle attribute arrays.

#include "cello.hpp"
#include "enzo.hpp"

EnzoSinkParticle::EnzoSinkParticle
(Block * block,
 int ib,
 int ip,
 int accretion_radius_cells,
 bool conserve_angular_momentum)
  : block_(block),
    batch_index_(ib),
    particle_index_(ip),
    accretion_radius_cells_(accretion_radius_cells),
    conserve_angular_momentum_(conserve_angular_momentum),
    total_mass_change_(0.0),
    total_momentum_x_change_(0.0),
    total_momentum_y_change_(0.0),
    total_momentum_z_change_(0.0),
    total_metal_mass_change_(0.0)
{
   // Read in the particle data
   Particle particle = block_->data()->particle();
   int it = particle.type_index("sink");
   bool metals = particle.has_attribute(it,"metal_fraction");
   enzo_float *pmass, *px, *py, *pz, *pvx, *pvy, *pvz, *paccrate, *pmetalfrac;

   const int ia_m       = particle.attribute_index (it, "mass");
   const int ia_x       = particle.attribute_index (it, "x");
   const int ia_y       = particle.attribute_index (it, "y");
   const int ia_z       = particle.attribute_index (it, "z");
   const int ia_vx      = particle.attribute_index (it, "vx");
   const int ia_vy      = particle.attribute_index (it, "vy");
   const int ia_vz      = particle.attribute_index (it, "vz");
   const int ia_accrate = particle.attribute_index (it, "accretion_rate");
   const int ia_mf      = particle.attribute_index (it, "metal_fraction");

   pmass      = (enzo_float *) particle.attribute_array(it, ia_m, batch_index_);
   px         = (enzo_float *) particle.attribute_array(it, ia_x, batch_index_);
   py         = (enzo_float *) particle.attribute_array(it, ia_y, batch_index_);
   pz         = (enzo_float *) particle.attribute_array(it, ia_z, batch_index_);
   pvx        = (enzo_float *) particle.attribute_array(it, ia_vx, batch_index_);
   pvy        = (enzo_float *) particle.attribute_array(it, ia_vy, batch_index_);
   pvz        = (enzo_float *) particle.attribute_array(it, ia_vz, batch_index_);
   paccrate   = (enzo_float *) particle.attribute_array(it, ia_accrate, batch_index_);
   pmetalfrac = metals ? (enzo_float *) particle.attribute_array(it, ia_mf, batch_index_)
                       : nullptr ;

   const int dm       = particle.stride(it, ia_m);
   const int dp       = particle.stride(it, ia_x);
   const int dv       = particle.stride(it, ia_vx);
   const int daccrate = particle.stride(it, ia_accrate);
   const int dmf      = metals ? particle.stride(it, ia_mf) : 0;

   mass_           = pmass[particle_index_ * dm];
   x_              = px[particle_index_ * dp];
   y_              = py[particle_index_ * dp];
   z_              = pz[particle_index_ * dp];
   vx_             = pvx[particle_index_ * dv];
   vy_             = pvy[particle_index_ * dv];
   vz_             = pvz[particle_index_ * dv];
   accretion_rate_ = paccrate[particle_index_ * daccrate];
   metal_fraction_ = metals ? pmetalfrac[particle_index_ * dmf] : 0.0;


   // Set accretion radius to be accretion_radius_cells_ multipled by minimum cell width
   double hx, hy, hz;
   block_->cell_width(&hx, &hy, &hz);
   const double min_cell_width = std::min(hx,std::min(hy,hz));
   accretion_radius_ = accretion_radius_cells_ * min_cell_width;

   // Find the bounding region of the accretion zone
   double xm, ym, zm;
   block_->data()->lower(&xm,&ym,&zm);
   int gx, gy, gz;
   block_->data()->field().ghost_depth(0,&gx,&gy,&gz);

   min_ind_x_ = ceil((x_ - xm - accretion_radius_) / hx - 0.5) + gx;
   min_ind_y_ = ceil((y_ - ym - accretion_radius_) / hy - 0.5) + gy;
   min_ind_z_ = ceil((z_ - zm - accretion_radius_) / hz - 0.5) + gz;
   max_ind_x_ = floor((x_ - xm + accretion_radius_) / hx - 0.5) + gx;
   max_ind_y_ = floor((y_ - ym + accretion_radius_) / hy - 0.5) + gy;
   max_ind_z_ = floor((z_ - zm + accretion_radius_) / hz - 0.5) + gz;
}

// -------------------------------------------------------------------------------------------

bool EnzoSinkParticle::cell_in_accretion_zone(int i, int j, int k,
					      double* norm_disp_x,
					      double* norm_disp_y,
					      double* norm_disp_z) throw() {


  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  double xm, ym, zm;
  block_->data()->lower(&xm,&ym,&zm);
  int gx, gy, gz;
  block_->data()->field().ghost_depth(0,&gx,&gy,&gz);

  // Find the coordinates of the center of the cell
  const double cell_center_x = xm + (i - gx) * hx;
  const double cell_center_y = ym + (j - gy) * hy;
  const double cell_center_z = zm + (k - gz) * hz;

  // Get the components of the displacement vector of the center of the cell from
  // the particle
  const double disp_x = cell_center_x - x_;
  const double disp_y = cell_center_y - y_;
  const double disp_z = cell_center_z - z_;

  // Compute the square of the magnitude of this vector
  const double r2 = disp_x * disp_x + disp_y * disp_y + disp_z * disp_z;

  // If angular momentum of gas is conserved, we need to compute the components
  // of the normalized displacement vector
  if (conserve_angular_momentum_){
     const double r = sqrt(r2);
     *norm_disp_x = disp_x / r;
     *norm_disp_y = disp_y / r;
     *norm_disp_z = disp_z / r;
  }

  // Return whether or not cell is in accretion zone
  return (r2 < accretion_radius_ * accretion_radius_);

}

// ---------------------------------------------------------------------------------------------

void EnzoSinkParticle::update_quantities(enzo_float density_change,
					 int index,
					 double norm_disp_x,
					 double norm_disp_y,
					 double norm_disp_z) throw() {

  int it = cello::particle_descr()->type_index("sink");

  // Get pointers to field data
  Field field = block_->data()->field();

  enzo_float * density       = (enzo_float*) field.values("density");
  enzo_float * vx            = (enzo_float*) field.values("velocity_x");
  enzo_float * vy            = (enzo_float*) field.values("velocity_y");
  enzo_float * vz            = (enzo_float*) field.values("velocity_z");
  enzo_float * metal_density =
    cello::particle_descr()->has_attribute(it,"metal_fraction") ?
    (enzo_float*) field.values("metal_density") : nullptr;

  enzo_float * density_source    = (enzo_float*) field.values("density_source");
  enzo_float * mom_dens_x_source = (enzo_float*) field.values("mom_dens_x_source");
  enzo_float * mom_dens_y_source = (enzo_float*) field.values("mom_dens_y_source");
  enzo_float * mom_dens_z_source = (enzo_float*) field.values("mom_dens_z_source");
  enzo_float * te_dens_source    = (enzo_float*) field.values("te_dens_source");

  // Get the cell volume
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  const double cell_volume = hx * hy * hz;

  // Get the mass change from this cell and update the total mass change
  const enzo_float mass_change = density_change * cell_volume;
  total_mass_change_ += mass_change;

  // Update total metal mass change if required
  if (cello::particle_descr()->has_attribute(it,"metal_fraction"))
    total_metal_mass_change_ = (density_change / density[index]) * metal_density[index];

  // Set density_sink equal to minus the density change
  density_source[index] = -density_change;

  // These variables store the change of the particle's momentum due to accretion
  // from this cell
  enzo_float momentum_x_change;
  enzo_float momentum_y_change;
  enzo_float momentum_z_change;
  if (conserve_angular_momentum_){

    // Get radial component of velocity
    // Technically, what I call the "radial velocity"
    // is really the radial velocity in the frame of reference
    // where the particle is at rest, plus the particle's
    // velocity.
    const double vx_radial = (vx[index] - vx_) * norm_disp_x + vx_;
    const double vy_radial = (vy[index] - vy_) * norm_disp_y + vy_;
    const double vz_radial = (vz[index] - vz_) * norm_disp_z + vz_;

    // Only the radial momentum of gas is changed (due to mass loss).
    // This momentum is added to particle
    momentum_x_change = mass_change * vx_radial;
    momentum_y_change = mass_change * vy_radial;
    momentum_z_change = mass_change * vz_radial;

  } else {

    // In this case, we don't change the velocity of the gas, so
    // momentum change of particle just depends on the mass removed
    // from the gas
    momentum_x_change = mass_change * vx[index];
    momentum_y_change = mass_change * vy[index];
    momentum_z_change = mass_change * vz[index];

  }

  // Update total particle momentum change
  total_momentum_x_change_ += momentum_x_change;
  total_momentum_y_change_ += momentum_y_change;
  total_momentum_z_change_ += momentum_z_change;

  // Set "mom_dens_source" fields to minus the particle's momentum change
  mom_dens_x_source[index] = -momentum_x_change;
  mom_dens_y_source[index] = -momentum_y_change;
  mom_dens_z_source[index] = -momentum_z_change;

  // Also need to compute change in specific kinetic energy
  const double old_specific_ke =
    0.5 * density[index] * ( vx[index] * vx[index] +
			     vy[index] * vy[index] +
			     vz[index] * vz[index] );

  const double new_vx =
    (density[index] * vx[index] - momentum_x_change) / (density[index] - density_change);
  const double new_vy =
    (density[index] * vy[index] - momentum_y_change) / (density[index] - density_change);
  const double new_vz =
    (density[index] * vz[index] - momentum_z_change) / (density[index] - density_change);

  const double new_specific_ke =
    0.5 * (density[index] - density_change) * ( new_vx * new_vx +
						new_vy * new_vy +
						new_vz * new_vz );

  te_dens_source[index] = new_specific_ke - old_specific_ke;
  return;
}

// -----------------------------------------------------------------------------------------

void EnzoSinkParticle::write_particle_data() throw() {


  Particle particle = block_->data()->particle();
  int it = particle.type_index("sink");
  bool metals = particle.has_attribute(it,"metal_fraction");

  // Get pointers to particle data
  enzo_float *pmass, *pvx, *pvy, *pvz, *paccrate, *pmetalfrac;

  const int ia_m       = particle.attribute_index (it, "mass");
  const int ia_vx      = particle.attribute_index (it, "vx");
  const int ia_vy      = particle.attribute_index (it, "vy");
  const int ia_vz      = particle.attribute_index (it, "vz");
  const int ia_accrate = particle.attribute_index (it, "accretion_rate");
  const int ia_mf      = particle.attribute_index (it, "metal_fraction");

  pmass      = (enzo_float *) particle.attribute_array(it, ia_m, batch_index_);
  pvx        = (enzo_float *) particle.attribute_array(it, ia_vx, batch_index_);
  pvy        = (enzo_float *) particle.attribute_array(it, ia_vy, batch_index_);
  pvz        = (enzo_float *) particle.attribute_array(it, ia_vz, batch_index_);
  paccrate   = (enzo_float *) particle.attribute_array(it, ia_accrate, batch_index_);
  pmetalfrac = metals ? (enzo_float *) particle.attribute_array(it, ia_mf, batch_index_)
                      : nullptr ;

  const int dm       = particle.stride(it, ia_m);
  const int dv       = particle.stride(it, ia_vx);
  const int daccrate = particle.stride(it, ia_accrate);
  const int dmf      = metals ? particle.stride(it, ia_mf) : 0;

  const enzo_float old_momentum_x = mass_ * vx_;
  const enzo_float old_momentum_y = mass_ * vy_;
  const enzo_float old_momentum_z = mass_ * vz_;

  const enzo_float old_metal_mass = mass_ * metal_fraction_;

  // Set new values for particle attributes
  pmass[particle_index_ * dm] = mass_ + total_mass_change_;

  pvx[particle_index_ * dv] =
    (old_momentum_x + total_momentum_x_change_) / pmass[particle_index_ * dm];
  pvy[particle_index_ * dv] =
    (old_momentum_y + total_momentum_y_change_) / pmass[particle_index_ * dm];
  pvz[particle_index_ * dv] =
    (old_momentum_z + total_momentum_z_change_) / pmass[particle_index_ * dm];

  paccrate[particle_index_ * daccrate] = total_mass_change_ / block_->dt();

  if (pmetalfrac) pmetalfrac[particle_index_ * dmf] =
      (old_metal_mass + total_metal_mass_change_) / pmass[particle_index_ * dm];

  return;
}
