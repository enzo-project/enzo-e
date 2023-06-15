// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodFBNetDeposit.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Wed Jun 7 16:14:38 PDT 2023
/// @brief    [\ref Enzo] Declaration of EnzoMethodFBNetDeposit
///           method for depositing Pop III remnants, identified using StarNet,
///           onto the mesh.

#ifndef ENZO_ENZO_METHOD_FBNET_DEPOSIT_HPP
#define ENZO_ENZO_METHOD_FBNET_DEPOSIT_HPP

class EnzoMethodFBNetDeposit : public Method {

  /// @class    EnzoMethodFBNetDeposit
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve heat equation
  /// using forward Euler method

public: // interface

  /// Create a new EnzoMethodFBNetDeposit object
  EnzoMethodFBNetDeposit();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodFBNetDeposit);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodFBNetDeposit (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block ) throw();

  virtual std::string name () throw () 
  { return "fbnet_deposit"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

protected: // methods

  void call_compute_done( EnzoBlock * block ) throw();

protected: // attributes

};

#endif /* ENZO_ENZO_METHOD_FBNET_DEPOSIT_HPP */
