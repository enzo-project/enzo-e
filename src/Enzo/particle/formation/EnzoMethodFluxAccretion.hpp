// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFluxAccretion.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       26 May 2022
/// @brief      Computes accretion rates according to "Flux Accretion" method, as
///             described in Section 5.3 of Bleuler, A & Teyssier, R 2004; MNRAS, 445, 4015-4036


#ifndef ENZO_ENZO_METHOD_FLUX_ACCRETION
#define ENZO_ENZO_METHOD_FLUX_ACCRETION

class EnzoMethodFluxAccretion : public EnzoMethodAccretion {

public:

  // Constructor
  EnzoMethodFluxAccretion(double accretion_radius_cells,
				double density_threshold,
				double max_mass_fraction);

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodFluxAccretion);

  // Charm++ PUP::able declarations
  EnzoMethodFluxAccretion (CkMigrateMessage *m)
   : EnzoMethodAccretion (m)
   {  }

  /// Charm++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply method
  virtual void compute ( Block * block) throw();

  /// Name
  virtual std::string name () throw()
   { return "accretion";}

private:

  // methods

  void compute_(Block * block);

};


#endif // ENZO_ENZO_METHOD_FLUX_ACCRETION
