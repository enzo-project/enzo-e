
#ifndef ENZO_ENZO_METHOD_VLCT_HPP
#define ENZO_ENZO_METHOD_VLCT_HPP
class EnzoMethodVlct : public Method {

  /// @class    EnzoMethodVlct
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate ENZO's VL + CT MHD method

public: // interface

  /// Creae a new EnzoMethodVlct object
  EnzoMethodVlct();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodVlct);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodVlct (CkMigrateMessage *m)
    : Method (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
public: // virtual methods

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "vlct"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
