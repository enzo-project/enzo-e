// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulence.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:57 UTC 2014
/// @brief    [\ref Enzo] Initial conditions for turbulence simulations

#ifndef ENZO_ENZO_INITIAL_TURBULENCE_HPP
#define ENZO_ENZO_INITIAL_TURBULENCE_HPP

class EnzoInitialTurbulence : public Initial {

  /// @class    EnzoInitialTurbulence
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for turbulence simulations

public: // interface

  /// CHARM++ constructor
  EnzoInitialTurbulence(double density_initial,
			double pressure_initial,
			double temperature_initial,
			double gamma) throw() 
    : Initial(),
      density_initial_(density_initial),
      pressure_initial_(pressure_initial),
      temperature_initial_(temperature_initial)
  { }
  
  /// Constructor
  EnzoInitialTurbulence(int cycle, double time,
			double density_initial,
			double pressure_initial,
			double temperature_initial,
			double gamma) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialTurbulence);

  /// CHARM++ migration constructor
  EnzoInitialTurbulence(CkMigrateMessage *m) 
    : Initial (m), 
      density_initial_(0),
      pressure_initial_(0),
      temperature_initial_(0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (
   Block * block,
   const FieldDescr * field_descr,
   const Hierarchy * hierarchy
   ) throw();


private: // attributes

  double density_initial_;
  double pressure_initial_;
  double temperature_initial_;
  double gamma_;

};

#endif /* ENZO_ENZO_INITIAL_TURBULENCE_HPP */

