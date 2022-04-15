// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionComputeBondiHoyle.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @author     John Regan (john.regan@mu.ie)
/// @date       13 April 2022
/// @brief      Computes accretion rates according to Bondi-Hoyle model.
///             See Krumholz+ 2004, ApJ, 611, 399 for details.



#ifndef ENZO_ENZO_METHOD_ACCRETION_COMPUTE_BONDI_HOYLE
#define ENZO_ENZO_METHOD_ACCRETION_COMPUTE_BONDI_HOYLE

class EnzoMethodAccretionComputeBondiHoyle : public EnzoMethodAccretionCompute {
  
public:
  
  // Constructor
  EnzoMethodAccretionComputeBondiHoyle(double accretion_radius_cells,
				       double density_threshold,
				       double max_mass_fraction,
				       bool conserve_angular_momentum);

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodAccretionComputeBondiHoyle);

  // Charm++ PUP::able declarations
  EnzoMethodAccretionComputeBondiHoyle (CkMigrateMessage *m)
   : EnzoMethodAccretionCompute (m)
   {  }

  /// Charm++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply method
  virtual void compute ( Block * block) throw();

  /// Name
  virtual std::string name () throw()
   { return "accretion_compute";}
  
// protected: // methods

//   void compute_(Block * block);
  

// private:

//   // attributes

//   // The maximum fraction of mass that can be accreted from a cell
//   // in one timestep
//   double max_mass_fraction_;

//   // methods
//   // Given the position of a particle, this function returns the
//   // index of the cell containing the particle
//   const std::vector<int> get_host_cell_index__
//   (const Block * block, enzo_float x, enzo_float y, enzo_float z)
//     throw();


//   // Given the position of a particle, this function returns a
//   // vector containing the indices of cells in the particle's accretion zone.
//   std::pair<const std::vector<int>,const std::vector<double>> get_accretion_zone_
//   (const Block * block, enzo_float x, enzo_float y, enzo_float z)
//     throw();

//   // Computes the Bondi-Hoyle radius and returns its value in cgs units.
//   // See Krumholz paper for more details)
//   const double compute_bondi_hoyle_radius_cgs_
//   (const Block * block, int host_cell_index, enzo_float mass,
//    enzo_float vx, enzo_float vy, enzo_float vz) throw();

//   // Compute normalized weights (i.e., the sum of the weights is 1), for
//   // cells in the accretion zone.
//   // See Krumholz paper for more details
//   const std::vector<double> compute_weights_
//   (const Block * block,
//    const std::vector<double>& r2,
//    const double r_BH_cgs) throw();

//   // Computes the mass removed from each cell in the accretion zone,
//   // setting values for density_accreted field in the process,
//   // and returns the mass removed from each cell as a vector
//   // pmass_cgs is the particle mass in cgs units
//   const std::vector<enzo_float> compute_mass_removed_
//   (Block * block, const std::vector<enzo_float>& weights,
//    const std::vector<int>& acc_zone_indices,
//    const int host_cell_index,
//    const double r_BH_cgs, const double pmass_cgs) throw();

//   // Returns a vector of length 3, containing the x, y, and z components
//   // of the total momentum which is removed from the gas, and must be
//   // added to the particle
//   const std::vector<enzo_float> compute_momentum_removed_
//   (const Block * block, const std::vector<enzo_float>& mass_removed) throw();

};


#endif /* EnzoMethodAccretionComputeBondiHoyle */
