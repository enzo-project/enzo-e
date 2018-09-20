// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodComovingExpansion.cpp
/// @author   Britton Smith (bds006@ucsd.edu)
/// @date     Wed May 24 12:25:56 PDT 2017
/// @brief    Implements comoving expansion class

#ifndef ENZO_ENZO_METHOD_COMOVING_EXPANSION_HPP
#define ENZO_ENZO_METHOD_COMOVING_EXPANSION_HPP

class EnzoMethodComovingExpansion : public Method {

  /// @class    EnzoMethodComovingExpansion
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo comoving expansion terms

public: // interface

  /// Create a new EnzoMethodComovingExpansion object
  EnzoMethodComovingExpansion(bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodComovingExpansion);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodComovingExpansion (CkMigrateMessage *m)
    : Method (m),
      comoving_coordinates_(false)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual methods
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "comoving_expansion"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

private: // attributes

  bool comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_COMOVING_EXPANSION_HPP */

