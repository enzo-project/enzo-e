// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Declaration of EnzoMethodHeat
///           forward Euler solver for the heat equation

#ifndef ENZO_ENZO_METHOD_HEAT_HPP
#define ENZO_ENZO_METHOD_HEAT_HPP

class EnzoMethodHeat : public Method {

  /// @class    EnzoMethodHeat
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve heat equation
  /// using forward Euler method

public: // interface

  /// Create a new EnzoMethodHeat object
  EnzoMethodHeat(const FieldDescr *, double alpha, double courant);

  EnzoMethodHeat() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodHeat);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodHeat (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "heat"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  template <class T>
  void compute_ (Block * block, T * Unew ) const throw();

protected: // attributes

  /// Thermal diffusivity
  double alpha_;

  /// Courant safety number
  double courant_;
};

#endif /* ENZO_ENZO_METHOD_HEAT_HPP */
