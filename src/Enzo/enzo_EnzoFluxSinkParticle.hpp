/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoFluxSinkParticle.hpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   26 May 2022
/// @brief  Implementation of EnzoFluxSinkParticle, a class which derives
///         from EnzoSinkParticle. Contains additional data and methods used for
///         computing accretion rates according to the method described in
///         Bleuler, A & Teyssier, R 2004; MNRAS, 445, 4015-4036

#ifndef ENZO_ENZO_FLUX_SINK_PARTICLE
#define ENZO_ENZO_FLUX_SINK_PARTICLE

class EnzoFluxSinkParticle : public EnzoSinkParticle{

  /// @class   EnzoSinkParticle
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Class derived from EnzoSinkParticle.
  ///          Contains additional data and methods used for
  ///          computing accretion rates according to the method described in
  ///          Bleuler, A & Teyssier, R 2004; MNRAS, 445, 4015-4036


public:

  /// Constructor
  /// Arguments: `block`                     - valid pointer to a block
  ///            `ib`                        - batch index
  ///            `ip`                        - index of particle within batch
  ///            `accretion_radius`          - radius of accretion zone
  ///            `density_threshold`         - density threshold set by the accretion method
  EnzoFluxSinkParticle(Block * block,
		       int ib,
		       int ip,
		       double accretion_radius,
		       double density_threshold);

  /// Computes the change in density for each cell in the accretion zone and updates
  /// properties of the sink particle
  /// Density change in a given cell is the accretion rate over cell volume multiplied by
  /// the timestep multiplied by the fraction of mass in the whole accretion zone which
  /// is in the cell (i.e., the mass in the cell divided by the mass in the accretion zone).
  /// Density change is limited to `max_mass_fraction` multiplied by the initial density,
  /// and furthermore, the new density cannot be lower than the density threshold.
  void compute(double max_mass_fraction) throw();

private:

  /// Methods

  /// This method does the following:
  /// - Finds which cells are in the accretion zone.
  /// - Computes the sum of densities in the accretion zone
  /// - Computes the total mass flux into the accretion zone divided by the cell volume
  ///   (Equation 35)
  /// - Computes the accretion rate divided by cell volume (Equation 36)
  void init_() throw();

  /// Computes and returns the mass flux divided by the cell volume through a given cell,
  /// specified by `ix`,`iy`,`iz`. See Equation 35.
  double get_mass_flux_over_cell_volume_(int ix,int iy,int iz) throw();

  /// Computes and sets accretion rate divided by cell volume according to Equation 36.
  void set_accretion_rate_over_cell_volume_
  (const double total_mass_flux_over_cell_volume,
   const double mean_density) throw();

  /// Attributes

  /// Vector containing the 1d indices of the cells in the accretion zone.
  /// Computed by compute_acc_zone_indices_();
  std::vector<int> acc_zone_1d_indices_;

  /// Value of the accretion rate. Computed by set_accretion_rate_over_cell_volume().
  /// This does not necessarily determine the total mass increase of the sink particle,
  /// since the actual density change of each cell will be limited by
  /// max_mass_fraction_ and density_threshold_
  double accretion_rate_over_cell_volume_;

  /// The sum of densities of cells in the accretion zone.
  double sum_of_densities_;

  /// The density threshold set by the accretion method
  double density_threshold_;

};

#endif // ENZO_ENZO_FLUX_SINK_PARTICLE
