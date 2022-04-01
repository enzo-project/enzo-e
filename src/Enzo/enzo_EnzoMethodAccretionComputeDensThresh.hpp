// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionCompute.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       10 March 2022
/// @brief      Implementation of EnzoMethodAccretionComputeDensThresh, a class
///             from EnzoMethodAccretionCompute.
///             This method reduces the gas density in the accretion zone around
///             a star particle to a value set by density_threshold_,
///             and adds mass and momentum lost by the gas to the star particle.

#ifndef ENZO_ENZO_METHOD_ACCRETION_COMPUTE_DENS_THRESH
#define ENZO_ENZO_METHOD_ACCRETION_COMPUTE_DENS_THRESH

class EnzoMethodAccretionComputeDensThresh : public EnzoMethodAccretionCompute {

public:

  // Constructor
  EnzoMethodAccretionComputeDensThresh(double accretion_radius_cells,
				       double density_threshold);

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodAccretionComputeDensThresh);

  /// Charm++ PUP::able migration constructor
  EnzoMethodAccretionComputeDensThresh (CkMigrateMessage *m)
    : EnzoMethodAccretionCompute (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "accretion_compute"; }

  /// Not sure if this is needed
  virtual std::string particle_type () throw()
  { return "star";}
  
private:
  // methods

  void compute_(Block * block);

  // attributes

  // The "density threshold" parameter: within the accretion zone around a star
  // particle, the density is reduced to the treshold, and the mass and momentum
  // lost by the gas is added to the star
  double density_threshold_;

};

#endif // ENZO_ENZO_METHOD_ACCRETION_COMPUTE_DENS_THRESH
