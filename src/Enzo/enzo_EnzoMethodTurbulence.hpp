// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 00:31:13 UTC 2014
/// @brief    [\ref Enzo] Implementation of Enzo TURBULENCE hydro method

#ifndef ENZO_ENZO_METHOD_TURBULENCE_HPP
#define ENZO_ENZO_METHOD_TURBULENCE_HPP

class EnzoMethodTurbulence : public Method {

  /// @class    EnzoMethodTurbulence
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's TURBULENCE hydro method

public: // interface

  /// Create a new EnzoMethodTurbulence object
  EnzoMethodTurbulence();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTurbulence);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodTurbulence (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( CommBlock * comm_block) throw();

  /// Compute maximum timestep for this method
  virtual double timestep ( CommBlock * comm_block) throw();

};

#endif /* ENZO_ENZO_METHOD_TURBULENCE_HPP */
