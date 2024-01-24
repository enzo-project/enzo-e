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
 double accretion_radius)
  : block_(block),
    batch_index_(ib),
    particle_index_(ip),
    accretion_radius_(accretion_radius),
    total_pmass_change_(0.0),
    total_momentum_x_change_(0.0),
    total_momentum_y_change_(0.0),
    total_momentum_z_change_(0.0),
    total_pmetal_mass_change_(0.0)
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

   pmass_           = pmass[particle_index_ * dm];
   px_              = px[particle_index_ * dp];
   py_              = py[particle_index_ * dp];
   pz_              = pz[particle_index_ * dp];
   pvx_             = pvx[particle_index_ * dv];
   pvy_             = pvy[particle_index_ * dv];
   pvz_             = pvz[particle_index_ * dv];
   accretion_rate_ = paccrate[particle_index_ * daccrate];
   pmetal_fraction_ = metals ? pmetalfrac[particle_index_ * dmf] : 0.0;

   // Find the bounding region of the accretion zone
   double xm, ym, zm;
   block_->data()->lower(&xm,&ym,&zm);
   int gx, gy, gz;
   block_->data()->field().ghost_depth(0,&gx,&gy,&gz);
   double hx, hy, hz;
   block->cell_width(&hx, &hy, &hz);

   min_ind_x_ = ceil((px_ - xm - accretion_radius_) / hx - 0.5) + gx;
   min_ind_y_ = ceil((py_ - ym - accretion_radius_) / hy - 0.5) + gy;
   min_ind_z_ = ceil((pz_ - zm - accretion_radius_) / hz - 0.5) + gz;
   max_ind_x_ = floor((px_ - xm + accretion_radius_) / hx - 0.5) + gx;
   max_ind_y_ = floor((py_ - ym + accretion_radius_) / hy - 0.5) + gy;
   max_ind_z_ = floor((pz_ - zm + accretion_radius_) / hz - 0.5) + gz;
}

// -------------------------------------------------------------------------------------------

bool EnzoSinkParticle::cell_in_accretion_zone(int i, int j, int k) throw()
{
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  double xm, ym, zm;
  block_->data()->lower(&xm,&ym,&zm);
  int gx, gy, gz;
  block_->data()->field().ghost_depth(0,&gx,&gy,&gz);

  // Find the coordinates of the center of the cell
  const double cell_center_x = xm + (i - gx + 0.5) * hx;
  const double cell_center_y = ym + (j - gy + 0.5) * hy;
  const double cell_center_z = zm + (k - gz + 0.5) * hz;

  // Get the components of the displacement vector of the center of the cell from
  // the particle
  const double disp_x = cell_center_x - px_;
  const double disp_y = cell_center_y - py_;
  const double disp_z = cell_center_z - pz_;

  // Compute the square of the magnitude of this vector
  const double r2 = disp_x * disp_x + disp_y * disp_y + disp_z * disp_z;

  // Return whether or not cell is in accretion zone
  return (r2 < accretion_radius_ * accretion_radius_);

}

// -------------------------------------------------------------------------------------------

bool EnzoSinkParticle::cell_in_accretion_zone(int i, int j, int k,
					      double* r2) throw()
{
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  double xm, ym, zm;
  block_->data()->lower(&xm,&ym,&zm);
  int gx, gy, gz;
  block_->data()->field().ghost_depth(0,&gx,&gy,&gz);

  // Find the coordinates of the center of the cell
  const double cell_center_x = xm + (i - gx + 0.5) * hx;
  const double cell_center_y = ym + (j - gy + 0.5) * hy;
  const double cell_center_z = zm + (k - gz + 0.5) * hz;

  // Get the components of the displacement vector of the center of the cell from
  // the particle
  const double disp_x = cell_center_x - px_;
  const double disp_y = cell_center_y - py_;
  const double disp_z = cell_center_z - pz_;

  // Compute the square of the magnitude of this vector
  *r2 = disp_x * disp_x + disp_y * disp_y + disp_z * disp_z;

  // Return whether or not cell is in accretion zone
  return (*r2 < accretion_radius_ * accretion_radius_);

}

// ---------------------------------------------------------------------------------------------

void EnzoSinkParticle::update(enzo_float density_change, int index) throw() {

  int it = cello::particle_descr()->type_index("sink");

  // Get pointers to field data
  Field field = block_->data()->field();

  enzo_float * density     = (enzo_float*) field.values("density");
  enzo_float * vx_gas      = (enzo_float*) field.values("velocity_x");
  enzo_float * vy_gas      = (enzo_float*) field.values("velocity_y");
  enzo_float * vz_gas      = (enzo_float*) field.values("velocity_z");

  enzo_float * metal_density =
    cello::particle_descr()->has_attribute(it,"metal_fraction") ?
    (enzo_float*) field.values("metal_density") : nullptr;

  enzo_float * density_source    = (enzo_float*) field.values("density_source");
  enzo_float * mom_dens_x_source = (enzo_float*) field.values("mom_dens_x_source");
  enzo_float * mom_dens_y_source = (enzo_float*) field.values("mom_dens_y_source");
  enzo_float * mom_dens_z_source = (enzo_float*) field.values("mom_dens_z_source");

  // Get the cell volume
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);
  const double cell_volume = hx * hy * hz;

  // Get the mass change from this cell and update the total mass change
  const enzo_float mass_change = density_change * cell_volume;
  total_pmass_change_ += mass_change;

  // Update total metal mass change if required
  if (cello::particle_descr()->has_attribute(it,"metal_fraction"))
    total_pmetal_mass_change_ = (density_change / density[index]) * metal_density[index];

  // Set density_sink equal to minus the density change
  density_source[index] = -density_change;

  // These variables store the change of the particle's momentum due to accretion
  // from this cell
  enzo_float momentum_x_change;
  enzo_float momentum_y_change;
  enzo_float momentum_z_change;

  // Compute change in momentum of particle due to accretion from this cell
  momentum_x_change = mass_change * vx_gas[index];
  momentum_y_change = mass_change * vy_gas[index];
  momentum_z_change = mass_change * vz_gas[index];

  // Update total particle momentum change
  total_momentum_x_change_ += momentum_x_change;
  total_momentum_y_change_ += momentum_y_change;
  total_momentum_z_change_ += momentum_z_change;

  // Set "mom_dens_source" fields to minus the particle's
  // momentum change divided by cell volume
  mom_dens_x_source[index] = -momentum_x_change / cell_volume;
  mom_dens_y_source[index] = -momentum_y_change / cell_volume;
  mom_dens_z_source[index] = -momentum_z_change / cell_volume;

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

  const enzo_float old_momentum_x = pmass_ * pvx_;
  const enzo_float old_momentum_y = pmass_ * pvy_;
  const enzo_float old_momentum_z = pmass_ * pvz_;

  const enzo_float old_metal_mass = pmass_ * pmetal_fraction_;

  // Set new values for particle attributes
  pmass[particle_index_ * dm] = pmass_ + total_pmass_change_;

  pvx[particle_index_ * dv] =
    (old_momentum_x + total_momentum_x_change_) / pmass[particle_index_ * dm];
  pvy[particle_index_ * dv] =
    (old_momentum_y + total_momentum_y_change_) / pmass[particle_index_ * dm];
  pvz[particle_index_ * dv] =
    (old_momentum_z + total_momentum_z_change_) / pmass[particle_index_ * dm];

  paccrate[particle_index_ * daccrate] = total_pmass_change_ / block_->state()->dt();

  if (pmetalfrac) pmetalfrac[particle_index_ * dmf] =
      (old_metal_mass + total_pmetal_mass_change_) / pmass[particle_index_ * dm];

  return;
}
