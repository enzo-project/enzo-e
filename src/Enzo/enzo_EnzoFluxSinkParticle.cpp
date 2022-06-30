/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoFluxSinkParticle.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   26 May 2022
/// @brief  Implementation of EnzoFluxSinkParticle, a class which derives
///         from EnzoSinkParticle. Contains additional data and methods used for
///         computing accretion rates according to the method described in
///         Bleuler, A & Teyssier, R 2004; MNRAS, 445, 4015-4036

#include "cello.hpp"
#include "enzo.hpp"

//-------------------------------------------------------------------------------

EnzoFluxSinkParticle::EnzoFluxSinkParticle
(Block * block,
 int ib,
 int ip,
 double accretion_radius,
 double density_threshold)
  : EnzoSinkParticle(block,ib,ip,accretion_radius),
    sum_of_densities_(0.0),
    density_threshold_(density_threshold)
{
  init_();
}

//--------------------------------------------------------------------------------

void EnzoFluxSinkParticle::init_() throw()
{
  // Get the field dimensions
  Field field = block_->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  // Get pointer to density field data
  enzo_float * density = (enzo_float*) field.values("density");

  // Will create a vector of indices of cells in the accretion zone.
  decltype(acc_zone_1d_indices_.size()) current_vector_index = 0;

  // Variable to store the total mass flux into the accretion zone, divided by the cell
  // volume.
  double total_mass_flux_over_cell_volume = 0.0;

  // Loop over cells in region bounding accretion zone
  for (int iz = min_ind_z_; iz < max_ind_z_; iz++){
    for (int iy = min_ind_y_; iy < max_ind_y_; iy++){
      for (int ix = min_ind_x_; ix < max_ind_x_; ix++){

	if (cell_in_accretion_zone(ix,iy,iz)) {
	  const int index = INDEX(ix,iy,iz,mx,my);

	  // append index to vector of accretion zone cell indices
	  acc_zone_1d_indices_.push_back(index);

	  // compute the mass flux over cell volume through this cell
	  // and increment total_mass_flux_over_cell_volume (Equation 35)
	  total_mass_flux_over_cell_volume += get_mass_flux_over_cell_volume_(ix,iy,iz);

	  // increment sum_of_densities_
	  sum_of_densities_ += density[index];

	} // if cell in accretion zone
      }
    }
  } // Triple loop over cells in region bounding accretion zone

  const double mean_density = sum_of_densities_ / acc_zone_1d_indices_.size();

  // Compute and set accretion rate over cell volume
  set_accretion_rate_over_cell_volume_(total_mass_flux_over_cell_volume, mean_density);

  return;
}

//--------------------------------------------------------------------------------------------

double EnzoFluxSinkParticle::get_mass_flux_over_cell_volume_(int ix,int iy,int iz) throw()
{
  // Get cell widths
  double hx, hy, hz;
  block_->cell_width(&hx, &hy, &hz);

  // Get field dimensions
  Field field = block_->data()->field();
  int mx, my, mz;
  field.dimensions(0, &mx, &my, &mz);

  // Get field data
  enzo_float * density = (enzo_float*) field.values("density");
  enzo_float * vx_gas  = (enzo_float*) field.values("velocity_x");
  enzo_float * vy_gas  = (enzo_float*) field.values("velocity_y");
  enzo_float * vz_gas  = (enzo_float*) field.values("velocity_z");

  // Compute divergence of "density * (gas velocity - particle velocity)" (Equation 36)
  const double d_rho_vx_dx =
    (density[INDEX(ix+1,iy,iz,mx,my)] * (vx_gas[INDEX(ix+1,iy,iz,mx,my)] - pvx_) -
     density[INDEX(ix-1,iy,iz,mx,my)] * (vx_gas[INDEX(ix-1,iy,iz,mx,my)] - pvx_)) / (2.0 * hx);
  const double d_rho_vy_dy =
    (density[INDEX(ix,iy+1,iz,mx,my)] * (vy_gas[INDEX(ix,iy+1,iz,mx,my)] - pvy_) -
     density[INDEX(ix,iy-1,iz,mx,my)] * (vy_gas[INDEX(ix,iy-1,iz,mx,my)] - pvy_)) / (2.0 * hy);
  const double d_rho_vz_dz =
    (density[INDEX(ix,iy,iz+1,mx,my)] * (vz_gas[INDEX(ix,iy,iz+1,mx,my)] - pvz_) -
     density[INDEX(ix,iy,iz-1,mx,my)] * (vz_gas[INDEX(ix,iy,iz-1,mx,my)] - pvz_)) / (2.0 * hz);

  return -(d_rho_vx_dx + d_rho_vy_dy + d_rho_vz_dz);

}

//--------------------------------------------------------------------------------------------

void EnzoFluxSinkParticle::set_accretion_rate_over_cell_volume_
(const double total_mass_flux_over_cell_volume,
 const double mean_density) throw()
{
  // Equation 36 - accretion rate is forced to be non-negative
  accretion_rate_over_cell_volume_ =
    std::max(0.0,(1.0 + 0.1*log(mean_density / density_threshold_))
	     * total_mass_flux_over_cell_volume);
  return;

}

//-------------------------------------------------------------------------------------

void EnzoFluxSinkParticle::compute(double max_mass_fraction) throw()
{
  // Check whether vectors have non-zero size and are all the same size
  ASSERT("EnzoFluxSinkParticle::compute_()",
	 "Function has been called when acc_zone_1d_indices_ has size of zero.",
	 acc_zone_1d_indices_.size() != 0);

  const int n_cells = acc_zone_1d_indices_.size();

  // Get pointer to density field data
  Field field = block_->data()->field();
  enzo_float * density = (enzo_float*) field.values("density");

  const double inv_sum_of_densities = 1.0 / sum_of_densities_;

  // Loop over cells in accretion zone
  for (int i = 0; i < n_cells; i++){
    const int index = acc_zone_1d_indices_[i];

    // Equation 37
    // Note: we have a density change on the LHS rather than a mass change, which is why
    // we have "accretion rate over cell volume" on the RHS rather than "accretion rate",
    // and instead of dividing by "n_cells times rho bar" we multiply by
    // "inverse sum of densities".
    double density_change =
      accretion_rate_over_cell_volume_ * density[index] * inv_sum_of_densities * block_->dt();

    // Restrict density change via density_threshold and max_mass_fraction
    density_change = std::min(std::min(density_change,density[index] - density_threshold_),
			      max_mass_fraction * density[index]);

    update(density_change,index);

  } // Loop over cells in accretion zone

  return;
}
