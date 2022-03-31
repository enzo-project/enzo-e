// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretionRemoveGas.hpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_ACCRETION_REMOVE_GAS
#define ENZO_ENZO_METHOD_ACCRETION_REMOVE_GAS

class EnzoMethodAccretionRemoveGas : public Method {

  /// @class   EnzoMethodAccretionRemoveGas
  /// @ingroup Enzo
  /// @brief   [\ref Enzo] Method for removing accreted gas from the gas density
  ///          field

public:

  /// Constructor
  EnzoMethodAccretionRemoveGas();

  /// Destructor
  virtual ~EnzoMethodAccretionRemoveGas() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodAccretionRemoveGas);

  /// Charm++ PUP::able migration constructor
  EnzoMethodAccretionRemoveGas (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method
  virtual void compute( Block * block) throw();

  /// Name
  virtual std::string name () throw()
  { return "accretion_remove_gas"; }

  // Compute the maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  void compute_(Block * block);

};

#endif /* EnzoMethodAccretionRemoveGas */
