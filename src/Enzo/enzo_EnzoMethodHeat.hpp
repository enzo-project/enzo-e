// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHeat.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Declaration of EnzoMethodForwardEuler
///           forward Euler

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
  EnzoMethodHeat(double alpha);

  EnzoMethodHeat() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodHeat);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodHeat (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute_block(FieldDescr *, CommBlock * block) throw();

  template <class T>
  void compute 
  (T * Unew,
   int ndx, int ndy, int ndz,
   int nx,  int ny,  int nz,
   double dt, double dx, double dy, double dz) const throw();

protected: // attributes

  double alpha_;

};

#endif /* ENZO_ENZO_METHOD_HEAT_HPP */
