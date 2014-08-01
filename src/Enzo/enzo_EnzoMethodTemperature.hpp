// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTemperature.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 15:53:28 PDT 2014
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeTemperature functions

#ifndef ENZO_ENZO_METHOD_TEMPERATURE_HPP
#define ENZO_ENZO_METHOD_TEMPERATURE_HPP

class EnzoMethodTemperature : public Method {

  /// @class    EnzoMethodTemperature
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputeTemperature functions

public: // interface

  /// Create a new EnzoMethodTemperature object
  EnzoMethodTemperature 
  (double density_floor,
   double temperature_floor,
   double mol_weight);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTemperature);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodTemperature (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
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

};

#endif /* ENZO_ENZO_METHOD_TEMPERATURE_HPP */
