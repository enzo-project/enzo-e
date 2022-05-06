/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretion.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   5 May 2022
/// @brief  Implementation of EnzoSinkParticle, a class which encapsulates
///         the data associated with a sink particle and its accretion zone,
///         as well as associated methods for computing accretion rates and
///         reading and writing to / from the particle attribute arrays.

#ifndef ENZO_ENZO_SINK_PARTICLE
#define ENZO_ENZO_SINK_PARTICLE

class EnzoSinkParticle {

  /// @class   EnzoSinkParticle
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Encapsulates data associated with a sink particle and its
  ///          accretion zone, as well as methods for computing accretion rates and
  ///          reading and writing to / from the particle attribute arrays.

public:

  /// Constructor
  /// Arguments: block                     - valid pointer to a block
  ///            ib                        - batch index
  ///            ip                        - index of particle within batch
  ///            accretion_radius_cells    - radius of accretion zone in terms of the cell
  ///                                      - width at the highest level of refinement.
  ///            conserve_angular_momentum - whether angular momentum of gas is conserved
  ///                                        during accretion
  EnzoSinkParticle(Block * block,
		   int ib,
		   int ip,
		   int accretion_radius_cells,
		   bool conserve_angular_momentum);

  /// Destructor
  virtual ~EnzoSinkParticle() throw() {};


  /// Returns whether cell with indices i,j,k is in the accretion zone.
  /// If angular momentum of gas is conserved, this function also 
  /// sets norm_disp_x/y/z to be the x/y/z components of the normalised
  /// displacement vector of the cell center relative to the particle
  /// position, i.e., it is a unit vector pointing from the particle to the
  /// center of the cell.
  bool cell_in_accretion_zone(int i, int j, int k,
			      double* norm_disp_x,
			      double* norm_disp_y,
			      double* norm_disp_z) throw();

  /// Takes as arguments the density change in a cell and its index, and the components
  /// of the normalised displacement vector of the cell center relative to the particle
  /// position, although if conserve_angular_momentum_ is false, these arguments are not
  /// used. This function updates the sink particle data, specifically, the total mass
  /// and momentum changes, and computes values for the source fields for the cell
  /// specified by `index`.
  void update_quantities(enzo_float density_change,
			 int index,
			 double norm_disp_x,
			 double norm_disp_y,
			 double norm_disp_z) throw();

  /// Writes particle data to the attribute arrays
  void write_particle_data() throw();

  /// Getter methods for the bounding indices
  inline int min_ind_x() {return min_ind_x_;}
  inline int min_ind_y() {return min_ind_y_;}
  inline int min_ind_z() {return min_ind_z_;}
  inline int max_ind_x() {return max_ind_x_;}
  inline int max_ind_y() {return max_ind_y_;}
  inline int max_ind_z() {return max_ind_z_;}

protected:

  /// Attributes

  /// Pointer to the block containing this particle
  Block * block_;

  /// The index of the batch containing this particle
  const int batch_index_;

  /// The index of the particle within the batch
  const int particle_index_;

  /// The radius of the accretion zone in terms of the cell width
  /// at the highest level of refinement
  const double accretion_radius_cells_;

  /// The accretion radius in code units.
  double accretion_radius_;

  /// Whether angular momentum is conserved
  const bool conserve_angular_momentum_;

  /// Cell indices which bound the accretion zone
  int min_ind_x_, max_ind_x_;
  int min_ind_y_, max_ind_y_;
  int min_ind_z_, max_ind_z_;

  /// Physical attributes of the particle
  enzo_float mass_;
  enzo_float x_, y_, z_;
  enzo_float vx_, vy_, vz_;
  enzo_float metal_fraction_;
  enzo_float accretion_rate_;

  /// Total changes of physical attributes of sink particle
  /// due to accretion
  enzo_float total_mass_change_;
  enzo_float total_momentum_x_change_;
  enzo_float total_momentum_y_change_;
  enzo_float total_momentum_z_change_;
  enzo_float total_metal_mass_change_;
};

#endif // ENZO_ENZO_SINK_PARTICLE
