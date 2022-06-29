/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoSinkParticle.hpp
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
  /// Arguments: `block`                     - valid pointer to a block
  ///            `ib`                        - batch index
  ///            `ip`                        - index of particle within batch
  ///            `accretion_radius`          - radius of accretion zone
  EnzoSinkParticle(Block * block,
		   int ib,
		   int ip,
		   double accretion_radius);

  /// Destructor
  virtual ~EnzoSinkParticle() throw() {};

  /// `i`,`j`,`k` are the (3D) indices of a cell.
  ///
  /// Function computes whether cell is inside accretion zone of sink particle.
  bool cell_in_accretion_zone(int i, int j, int k) throw();

  /// `i`,`j`,`k` are the (3D) indices of a cell.
  ///
  /// Function computes whether cell is inside accretion zone of sink particle.
  /// This version also sets r2 to be equal to the square of the distance of the
  /// center of the cell from the sink particle.
  bool cell_in_accretion_zone(int i, int j, int k, double* r2) throw();

  /// `density_change` is the change in density in given cell (specified by `index`)
  /// due to accretion.
  ///
  /// `index` is the (1D) index of the cell.
  ///
  /// This function updates the sink particle date and computes values for the source fields
  /// in the given cell (specified by `index`).
  void update(enzo_float density_change, int index) throw();

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

  /// The radius of the accretion zone.
  const double accretion_radius_;

  /// Cell indices which bound the accretion zone
  int min_ind_x_, max_ind_x_;
  int min_ind_y_, max_ind_y_;
  int min_ind_z_, max_ind_z_;

  /// Physical attributes of the particle
  enzo_float pmass_;
  enzo_float px_, py_, pz_;
  enzo_float pvx_, pvy_, pvz_;
  enzo_float pmetal_fraction_;
  enzo_float accretion_rate_;

  /// Total changes of physical attributes of sink particle
  /// due to accretion
  enzo_float total_pmass_change_;
  enzo_float total_momentum_x_change_;
  enzo_float total_momentum_y_change_;
  enzo_float total_momentum_z_change_;
  enzo_float total_pmetal_mass_change_;
};

#endif // ENZO_ENZO_SINK_PARTICLE
