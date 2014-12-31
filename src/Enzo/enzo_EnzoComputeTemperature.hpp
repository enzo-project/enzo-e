// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeTemperature.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeTemperature functions

#ifndef ENZO_ENZO_COMPUTE_TEMPERATURE_HPP
#define ENZO_ENZO_COMPUTE_TEMPERATURE_HPP

class EnzoComputeTemperature : public Compute {

  /// @class    EnzoComputeTemperature
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputeTemperature functions

public: // interface

  /// Create a new EnzoComputeTemperature object
  EnzoComputeTemperature 
  (double density_floor,
   double temperature_floor,
   double mol_weight,
   int comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeTemperature);
  
  /// Charm++ PUP::able migration constructor
  EnzoComputeTemperature (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Perform the computation on the block
  virtual void compute( CommBlock * comm_block) throw();

private: // functions

  template <typename T>
  void compute_(CommBlock * comm_block);

private: // attributes

  // minimum density
  double density_floor_;

  // minimum temperature: default 1.0
  double temperature_floor_;

  // mol weight: default 0.6
  double mol_weight_;

  int comoving_coordinates_;

};

#endif /* ENZO_ENZO_COMPUTE_TEMPERATURE_HPP */
