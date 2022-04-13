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
  EnzoMethodAccretionComputeBondiHoyle(double accretion_radius_cells);

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
  
protected: // methods

  void compute_(Block * block);
  
};

#endif /* EnzoMethodAccretionComputeBondiHoyle */
