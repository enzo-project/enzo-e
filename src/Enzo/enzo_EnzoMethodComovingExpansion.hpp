// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodComovingExpansion.cpp
/// @author   Britton Smith (bds006@ucsd.edu)
/// @date     Wed May 24 12:25:56 PDT 2017
/// @brief    Implements comoving expansion class

#ifndef ENZO_ENZO_METHOD_COMOVINGEXPANSION_HPP
#define ENZO_ENZO_METHOD_COMOVINGEXPANSION_HPP

class EnzoMethodComovingExpansion : public Method {

  /// @class    EnzoMethodComovingExpansion
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo comoving expansion terms

public: // interface

  /// Create a new EnzoMethodComovingExpansion object
  EnzoMethodComovingExpansion(const FieldDescr * field_descr,
                              EnzoConfig * enzo_config);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodComovingExpansion);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodComovingExpansion (CkMigrateMessage *m)
    : comoving_coordinates_(false)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "comoving_expansion"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

private: // attributes

  int comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_COMOVINGEXPANSION_HPP */

