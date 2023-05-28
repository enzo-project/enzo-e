// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodThresholdAccretion.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       10 March 2022
/// @brief      Implementation of EnzoMethodThresholdAccretion, a class
///             from EnzoMethodAccretion.
///             This method reduces the gas density in the accretion zone around
///             a sink particle to a value set by density_threshold_,
///             and adds mass and momentum lost by the gas to the sink particle.

#ifndef ENZO_ENZO_METHOD_THRESHOLD_ACCRETION
#define ENZO_ENZO_METHOD_THRESHOLD_ACCRETION

class EnzoMethodThresholdAccretion : public EnzoMethodAccretion {

public:

  // Constructor
  EnzoMethodThresholdAccretion(double accretion_radius_cells,
				double density_threshold,
				double max_mass_fraction);

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodThresholdAccretion);

  /// Charm++ PUP::able migration constructor
  EnzoMethodThresholdAccretion (CkMigrateMessage *m)
    : EnzoMethodAccretion (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "accretion"; }
  
private:
  // methods

  void compute_(Block * block);

};

#endif // ENZO_ENZO_METHOD_THRESHOLD_ACCRETION
