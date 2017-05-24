// See LICENSE_CELLO file for license and copyright information

// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#ifndef ENZO_ENZO_METHOD_COMOVINGEXPANSION_HPP
#define ENZO_ENZO_METHOD_COMOVINGEXPANSION_HPP

class Classname {

  /// @class    Classname
  /// @ingroup  Component
  /// @brief    [\ref Component] 

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

