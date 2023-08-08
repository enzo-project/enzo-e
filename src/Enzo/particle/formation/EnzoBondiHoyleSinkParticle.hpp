/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoBondiHoyleSinkParticle.hpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @author John Regan (john.regan@mu.ie)
/// @date   17 May 2022
/// @brief  Implementation of EnzoBondiHoyleSinkParticle, a class which derives
///         from EnzoSinkParticle. Contains additional data and methods used for
///         computing accretion rates according to the method described in
///         Krumholz+ 2004, ApJ, 611, 399

#ifndef ENZO_ENZO_BONDI_HOYLE_SINK_PARTICLE
#define ENZO_ENZO_BONDI_HOYLE_SINK_PARTICLE

class EnzoBondiHoyleSinkParticle : public EnzoSinkParticle{

  /// @class   EnzoSinkParticle
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Class derived from EnzoSinkParticle.
  ///          Contains additional data and methods used for
  ///          computing accretion rates according to the method described in
  ///          Krumholz+ 2004, ApJ, 611, 399


public:

  /// Constructor
  /// Arguments: `block`                     - valid pointer to a block
  ///            `ib`                        - batch index
  ///            `ip`                        - index of particle within batch
  ///            `accretion_radius`          - radius of accretion zone
  ///            `const_G`                   - value of gravitation constant in code units

  EnzoBondiHoyleSinkParticle(Block * block,
			     int ib,
			     int ip,
			     double accretion_radius,
			     double const_G);

  /// Computes the change in density for each cell in the accretion zone and updates
  /// properties of the sink particle
  /// Density change in a given cell is the cell weight multiplied by the bound fraction
  /// multiplied by the accretion rate multiplied by the timestep divided by the cell volume.
  /// Density change is limited to `max_mass_fraction` multiplied by the initial density,
  /// and furthermore, the new density cannot be lower than the density threshold.
  void compute(double density_threshold,
	       double max_mass_fraction) throw();

private:

  /// Methods

  /// This method does the following:
  /// - Finds the cell which contains the sink particle.
  /// - Computes the Bondi-Hoyle radius (Equation 10)
  /// - Finds which cells are in the accretion zone.
  /// - Computes the "bound fraction" for each cell (equation 15 of Krumholz paper and
  ///    accompanying text).
  /// - Computes the cell weights (Equations 13 and 14)
  /// - Computes the accretion rate (Equations 11 and 12)
  void init_() throw();

  /// Computes and sets the (1d and 3d) indices of the cell containing the sink particle
  void set_host_cell_indices_() throw();

  /// Computes and sets the Bondi-Hoyle radius according to Equation 10 of Krumholz paper
  void set_bondi_hoyle_radius_() throw();

  /// Computes and sets `v_inf_2_`
  void set_v_inf_2_() throw();

  /// Computes and sets `c_s_inf_2_`
  void set_c_s_inf_2_() throw();

  /// Computes and returns cell weight, given the square of the distance of the cell
  /// from the sink particle, according to Equations 13 and 14 of Krumholz paper
  double get_cell_weight_(double r2) throw();

  /// Computes and returns the bound fraction of a given cell, specified by `ix`,`iy`,`iz`,
  /// according Equation 15 of Krumholz paper and accompanying text.
  double get_cell_bound_fraction_(int ix, int iy, int iz) throw();

  /// Given coordinates `x`, `y`, `z`, computes values for `vx_rel`,`vy_rel`, and `vz_rel`
  /// by linear interpolation from the gas velocity field
  void compute_relative_velocity_
  (const double x,const double y,const double z,
   double* vx_rel, double* vy_rel,double* vz_rel) throw();

  /// Computes and sets accretion rate according to Equations 11 and 12 of Krumholz paper.
  void set_accretion_rate_() throw();

  /// "alpha function" that appears in Equation 12 of Krumholz paper.
  /// It is the dimensionless density profile which is a solution to simple Bondi accretion.
  double alpha_(double x) throw();

  /// Attributes

  /// The value of the gravitational constant in code units
  const double const_G_;

  /// 1d index of the cell containing the sink particle. Computed by compute_host_cell_index_();
  int host_cell_1d_index_;

  /// 3d indicies of the cell containing the sink;
  int host_cell_ix_, host_cell_iy_, host_cell_iz_;

  /// Vector containing the 1d indices of the cells in the accretion zone.
  /// Computed by compute_acc_zone_indices_();
  std::vector<int> acc_zone_1d_indices_;

  /// Vector of normalised weights (i.e., the sum of the weights is 1), for each cell
  /// in accretion zone. Computed by compute_cell_weight_()
  std::vector<double> cell_weights_;

  /// Vector of "bound fractions" for each cell. Used to reduce the fraction of mass
  /// accreted from each cell. Computed by compute_bound_fraction_();
  std::vector<double> bound_fractions_;

  /// Value of the Bondi-Hoyle radius. Computed by compute_bondi_hoyle_radius_()
  double r_BH_;

  /// The square of the magnitude of the velocity of gas in host cell relative to the sink
  /// particle.
  double v_inf_2_;

  /// The square of the sound speed of the gas in the host cell
  double c_s_inf_2_;

  /// Value of the accretion rate. Computed by set_accretion_rate_().
  /// This does not necessarily determine the total mass increase of the sink particle,
  /// since the actual density change of each cell will be limited by
  /// max_mass_fraction_ and density_threshold_
  double accretion_rate_;

};

#endif // ENZO_ENZO_BONDI_HOYLE_SINK_PARTICLE
