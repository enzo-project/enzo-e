// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodBondiHoyleAccretion.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @author     John Regan (john.regan@mu.ie)
/// @date       13 April 2022
/// @brief      Computes accretion rates according to Bondi-Hoyle model.
///             See Krumholz+ 2004, ApJ, 611, 399 for details.


#ifndef ENZO_ENZO_METHOD_BONDI_HOYLE_ACCRETION
#define ENZO_ENZO_METHOD_BONDI_HOYLE_ACCRETION

class EnzoMethodBondiHoyleAccretion : public EnzoMethodAccretion {
  
public:
  
  // Constructor
  EnzoMethodBondiHoyleAccretion(double accretion_radius_cells,
				double density_threshold,
				double max_mass_fraction);

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodBondiHoyleAccretion);

  // Charm++ PUP::able declarations
  EnzoMethodBondiHoyleAccretion (CkMigrateMessage *m)
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


#endif // ENZO_ENZO_METHOD_BONDI_HOYLE_ACCRETION
